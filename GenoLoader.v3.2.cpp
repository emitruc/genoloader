// =============================================================================
//  GenoLoader v3.2 — C++ translation (GCC/Linux)
//  Translated from GiNOLOADER_v3_2.py
//
//  Compile:
//      g++ -O2 -std=c++17 -o genoloader genoloader.cpp
//
//  Usage:
//      ./genoloader <vcf_file> [options]
//
//  Options:
//      --p1    <file>      Text file listing Pop1 sample names (default: none)
//      --p2    <file>      Text file listing Pop2 sample names (default: none)
//      --p0    <file>      Text file listing outgroup sample names (default: none)
//      --m1    <int>       Min individuals required in Pop1 (default: 0)
//      --m2    <int>       Min individuals required in Pop2 (default: 0)
//      --m0    <int>       Min individuals required in outgroup (default: 0)
//      --polX  <string>    Polarization strategy: REF | POP_OUT | POP1 | POP2 | POP1POP2 (default: REF)
//      --low_cov <YES|NO>  Enable low-coverage read resampling (default: NO)
//
//  Example:
//      ./genoloader sample.ann.vcf --p1 pop1.txt --p2 pop2.txt --p0 outgroup.txt
//                   --m1 5 --m2 5 --m0 2 --polX POP_OUT --low_cov YES
// =============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <map>

// =============================================================================
//  Utility: trim trailing whitespace / carriage returns from a string
// =============================================================================
static std::string rstrip(const std::string& s) {
    size_t end = s.find_last_not_of(" \t\r\n");
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// =============================================================================
//  Utility: split a string by a single-character delimiter
// =============================================================================
static std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> tokens;
    std::istringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delim))
        tokens.push_back(token);
    return tokens;
}

// =============================================================================
//  Utility: replace all occurrences of 'from' with 'to' in string 's'
// =============================================================================
static std::string replaceAll(std::string s, const std::string& from, const std::string& to) {
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
    return s;
}

// =============================================================================
//  Utility: check if string contains substring
// =============================================================================
static bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

// =============================================================================
//  read_file_2_list
//  Reads a plain-text file and returns non-empty trimmed lines as a vector.
// =============================================================================
static std::vector<std::string> read_file_2_list(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Cannot open file: " + filename);
    std::vector<std::string> result;
    std::string line;
    while (std::getline(f, line)) {
        std::string stripped = rstrip(line);
        if (!stripped.empty())
            result.push_back(stripped);
    }
    return result;
}

// =============================================================================
//  replace_geno
//  Inverts a genotype string: swaps REF (0) and ALT (1) alleles.
//  Used when the ancestral allele is the VCF ALT allele.
// =============================================================================
static std::string replace_geno(const std::string& text) {
    static const std::unordered_map<std::string, std::string> table = {
        {"0/0", "1|1"}, {"0/1", "1|0"}, {"1/0", "0|1"}, {"1/1", "0|0"},
        {"./0", "1|1"}, {"0/.", "1|1"}, {"./1", "0|0"}, {"1/.", "0|0"},
        {"1",   "0|0"}, {"0",   "1|1"}
    };
    auto it = table.find(text);
    return (it != table.end()) ? it->second : text;
}

// =============================================================================
//  replace_geno_coding
//  Converts a genotype string to derived-allele dosage: 0, 1, 2, or nan.
// =============================================================================
static std::string replace_geno_coding(const std::string& text) {
    static const std::unordered_map<std::string, std::string> table = {
        {"0/0", "0"}, {"0/1", "1"}, {"1/0", "1"}, {"1/1", "2"},
        {"./.", "nan"}, {".",   "nan"},
        {"./0", "0"}, {"0/.", "0"}, {"./1", "2"}, {"1/.", "2"},
        {"1",   "2"}, {"0",   "0"}
    };
    auto it = table.find(text);
    return (it != table.end()) ? it->second : text;
}

// =============================================================================
//  get_alleles
//  Splits each genotype string on '/' and returns a flat list of alleles.
//  Haploid genotypes (no '/') are duplicated.
// =============================================================================
static std::vector<std::string> get_alleles(const std::vector<std::string>& geno_samples) {
    std::vector<std::string> alleles;
    for (const auto& sample : geno_samples) {
        auto parts = split(sample, '/');
        if (parts.size() == 1) {
            alleles.push_back(parts[0]);
            alleles.push_back(parts[0]);
        } else {
            alleles.push_back(parts[0]);
            alleles.push_back(parts[1]);
        }
    }
    return alleles;
}

// =============================================================================
//  get_alleles_from_index
//  Extracts alleles for a given population from a VCF line using the header
//  index to locate each sample column.
// =============================================================================
static std::vector<std::string> get_alleles_from_index(
    const std::vector<std::string>& header,
    const std::vector<std::string>& line,
    const std::vector<std::string>& pop)
{
    std::vector<std::string> geno;
    for (const auto& sample : pop) {
        auto it = std::find(header.begin(), header.end(), sample);
        if (it == header.end())
            throw std::runtime_error("Sample not found in VCF header: " + sample);
        int idx = (int)(it - header.begin());
        // Extract GT field (first ':'-delimited token)
        auto fields = split(line[idx], ':');
        geno.push_back(fields[0]);
    }
    return get_alleles(geno);
}

// =============================================================================
//  mode_int
//  Returns the modal (most frequent) value from a vector of integers.
//  Ties are broken in favour of the smallest value (deterministic, matching
//  scipy.stats.mode behaviour: always returns 0 in a 50/50 split).
// =============================================================================
static int mode_int(const std::vector<int>& v) {
    std::map<int, int> counts;
    for (int x : v) counts[x]++;
    int best_val = v[0], best_cnt = 0;
    for (auto& kv : counts) {
        if (kv.second > best_cnt) {
            best_cnt = kv.second;
            best_val = kv.first;
        }
    }
    return best_val;
}

// =============================================================================
//  unique_count
//  Returns the number of distinct values in a vector.
// =============================================================================
static size_t unique_count(const std::vector<int>& v) {
    return std::set<int>(v.begin(), v.end()).size();
}

// =============================================================================
//  get_ancestral_allele
//  Core polarization logic. Given allele vectors for the outgroup (geno_out),
//  population 1 (geno1), and population 2 (geno2) — all as strings, with '.'
//  for missing — determines the ancestral allele (0 or 1) and assigns a
//  diagnostic flag describing the polarization scenario.
//
//  Returns: {ref_allele, flag}
// =============================================================================
static std::pair<int, std::string> get_ancestral_allele(
    const std::vector<std::string>& geno_out,
    const std::vector<std::string>& geno1,
    const std::vector<std::string>& geno2)
{
    // Filter missing values and convert to int
    auto to_int_vec = [](const std::vector<std::string>& g) {
        std::vector<int> v;
        for (const auto& x : g)
            if (x != ".") v.push_back(std::stoi(x));
        return v;
    };

    std::vector<int> out_geno = to_int_vec(geno_out);
    std::vector<int> p1_geno  = to_int_vec(geno1);
    std::vector<int> p2_geno  = to_int_vec(geno2);

    std::vector<int> p12_geno = p1_geno;
    p12_geno.insert(p12_geno.end(), p2_geno.begin(), p2_geno.end());

    std::vector<int> total_alleles = p12_geno;
    total_alleles.insert(total_alleles.end(), out_geno.begin(), out_geno.end());

    size_t n_out  = unique_count(out_geno);
    size_t n_p12  = unique_count(p12_geno);
    size_t n_p1   = unique_count(p1_geno);
    size_t n_p2   = unique_count(p2_geno);

    int ref_allele = 9;
    std::string flag;

    // --- All missing ---
    if (out_geno.empty() && p12_geno.empty()) {
        ref_allele = 9;
        flag = "allMiss";

    // --- Ingroup populations missing ---
    } else if (p12_geno.empty()) {
        ref_allele = mode_int(out_geno);
        flag = "inMiss";

    // --- Outgroup missing ---
    } else if (out_geno.empty() && !p12_geno.empty()) {

        if (p1_geno.empty()) {
            // pop1 missing
            ref_allele = mode_int(p2_geno);
            flag = (n_p2 == 1) ? "in2FoldAnc" : "in2FoldSeg";

        } else if (p2_geno.empty()) {
            // pop2 missing
            ref_allele = mode_int(p1_geno);
            flag = (n_p1 == 1) ? "in1FoldAnc" : "in1FoldSeg";

        } else if (n_p1 == 1 && n_p2 == 1) {
            ref_allele = p1_geno[0]; // conservative: use pop1
            flag = (p1_geno[0] == p2_geno[0]) ? "InvarOutMiss" : "altFix";

        } else if (n_p1 == 2 && n_p2 == 1) {
            ref_allele = p2_geno[0];
            flag = "unfoldFix2outMiss";

        } else if (n_p1 == 1 && n_p2 == 2) {
            ref_allele = p1_geno[0];
            flag = "unfoldFix1outMiss";

        } else if (n_p1 == 2 && n_p2 == 2) {
            ref_allele = mode_int(p12_geno); // not random; 0 wins ties
            flag = "inFold";
        }

    // --- All fixed (outgroup and ingroup monomorphic) ---
    } else if (n_out == 1 && n_p12 == 1) {
        ref_allele = p12_geno[0]; // conservative: use ingroup
        flag = "allFix";

    // --- Outgroup monomorphic, ingroup polymorphic (standard unfolded case) ---
    } else if (n_out == 1 && n_p12 == 2) {
        ref_allele = out_geno[0];
        if (p1_geno.empty() && !p2_geno.empty()) {
            flag = "unfolded2Miss1";
        } else if (!p1_geno.empty() && p2_geno.empty()) {
            flag = "unfolded1Miss2";
        } else if (n_p1 == 1 && n_p2 == 1) {
            flag = "unfoldedAltFix";
        } else if (n_p1 == 2 && n_p2 == 2) {
            flag = "unfoldedILS";
        } else if (n_p1 == 1 && n_p2 == 2) {
            flag = (p1_geno[0] == out_geno[0]) ? "unfoldedFixAnc1" : "unfoldedFixDer1";
        } else if (n_p1 == 2 && n_p2 == 1) {
            flag = (p2_geno[0] == out_geno[0]) ? "unfoldedFixAnc2" : "unfoldedFixDer2";
        }

    // --- Outgroup polymorphic, ingroup monomorphic ---
    } else if (n_out == 2 && n_p12 == 1) {
        ref_allele = p12_geno[0];
        flag = "InFixAnc";

    // --- Both outgroup and ingroup polymorphic (folded fallback) ---
    } else if (n_out == 2 && n_p12 == 2) {
        ref_allele = mode_int(total_alleles); // not random; 0 wins ties
        flag = "allFold";

    // --- Non-biallelic or unhandled case ---
    } else {
        ref_allele = 9;
        std::cerr << "WARNING: Non biallelic locus in the vcf. "
                     "Replaced with missing value for all samples\n";
        flag = "triall";
    }

    return {ref_allele, flag};
}

// =============================================================================
//  random_allele_sampling
//  For low-coverage data: samples one random read per individual from the
//  allelic depth (AD) field. Constructs a pool of alleles weighted by read
//  depth, shuffles it, and returns the first element.
//  Ancestral allele identity is used to map REF/ALT read counts to 0/1.
// =============================================================================
static std::vector<std::string> random_allele_sampling(
    const std::vector<std::string>& all_samples_coverage,
    int ancestral_allele,
    std::mt19937& rng)
{
    std::vector<std::string> allele_random;
    for (const auto& AD : all_samples_coverage) {
        if (AD == "." || AD == "0:0") {
            allele_random.push_back(".");
            continue;
        }
        auto parts = split(AD, ':');
        int ad0 = std::stoi(parts[0]);
        int ad1 = std::stoi(parts[1]);

        std::vector<std::string> pool;
        if (ancestral_allele == 0) {
            for (int i = 0; i < ad0; ++i) pool.push_back("0");
            for (int j = 0; j < ad1; ++j) pool.push_back("1");
        } else {
            for (int i = 0; i < ad0; ++i) pool.push_back("1");
            for (int j = 0; j < ad1; ++j) pool.push_back("0");
        }
        std::shuffle(pool.begin(), pool.end(), rng);
        allele_random.push_back(pool[0]);
    }
    return allele_random;
}

// =============================================================================
//  write_gt
//  Main function. Parses the annotated VCF, polarizes each SNP, recodes
//  genotypes as derived-allele dosage (0/1/2/nan), and writes output files.
// =============================================================================
static void write_gt(
    const std::string& vcf_file,
    const std::string& p1      = "",
    const std::string& p2      = "",
    const std::string& p0      = "",
    int m1 = 0, int m2 = 0, int m0 = 0,
    const std::string& polX    = "REF",
    const std::string& low_cov = "NO")
{
    // --- Load population sample lists ---
    std::vector<std::string> pop1, pop2, pop_out;

    if (p1.empty()) {
        std::cout << "Pop1 is empty\n";
    } else {
        pop1 = read_file_2_list(p1);
        std::cout << "Pop1 includes the following individuals:\n";
        for (size_t i = 0; i < pop1.size(); ++i)
            std::cout << pop1[i] << (i + 1 < pop1.size() ? "\t" : "\n");
    }
    if (p2.empty()) {
        std::cout << "Pop2 is empty\n";
    } else {
        pop2 = read_file_2_list(p2);
        std::cout << "Pop2 includes the following individuals:\n";
        for (size_t i = 0; i < pop2.size(); ++i)
            std::cout << pop2[i] << (i + 1 < pop2.size() ? "\t" : "\n");
    }
    if (p0.empty()) {
        std::cout << "Pop_out is empty\n";
    } else {
        pop_out = read_file_2_list(p0);
        std::cout << "Pop_out includes the following individuals:\n";
        for (size_t i = 0; i < pop_out.size(); ++i)
            std::cout << pop_out[i] << (i + 1 < pop_out.size() ? "\t" : "\n");
    }

    std::cout << "Minimum required number of individuals in pop1: "    << m1   << "\n";
    std::cout << "Minimum required number of individuals in pop2: "    << m2   << "\n";
    std::cout << "Minimum required number of individuals in pop_out: " << m0   << "\n";
    std::cout << "Polarization strategy is: "                          << polX << "\n";

    // --- Open output files ---
    std::string gt_filename = vcf_file;
    {
        size_t pos = gt_filename.rfind(".vcf");
        if (pos != std::string::npos)
            gt_filename.replace(pos, 4, "." + polX + ".gt");
        else
            gt_filename += "." + polX + ".gt";
    }
    std::ofstream gt(gt_filename);
    if (!gt.is_open())
        throw std::runtime_error("Cannot open output file: " + gt_filename);

    std::ofstream gt_cov;
    if (low_cov == "YES") {
        std::cout << "One random read resampling per individual genotype is switched on\n";
        std::string cov_filename = vcf_file;
        size_t pos = cov_filename.rfind(".vcf");
        if (pos != std::string::npos)
            cov_filename.replace(pos, 4, ".gt.covRandom1");
        else
            cov_filename += ".gt.covRandom1";
        gt_cov.open(cov_filename);
        if (!gt_cov.is_open())
            throw std::runtime_error("Cannot open low-cov output file: " + cov_filename);
    }

    // --- Counters ---
    long valid_loci   = 0;
    long inverted_loci = 0;
    long missing_loci = 0;

    // --- RNG for low-coverage resampling ---
    std::mt19937 rng(std::random_device{}());

    // --- Progress counter (replaces tqdm) ---
    long line_count = 0;
    {
        std::ifstream tmp(vcf_file);
        std::string dummy;
        while (std::getline(tmp, dummy)) ++line_count;
    }
    std::cout << "Parsing " << line_count << " lines...\n";

    // --- Parse VCF ---
    std::ifstream handle(vcf_file);
    if (!handle.is_open())
        throw std::runtime_error("Cannot open VCF file: " + vcf_file);

    std::vector<std::string> header;
    std::string raw_line;
    long processed = 0;
    int report_every = std::max(1L, line_count / 20); // report progress ~every 5%

    while (std::getline(handle, raw_line)) {
        ++processed;
        if (processed % report_every == 0)
            std::cout << "  " << (100 * processed / line_count) << "% done\r" << std::flush;

        std::string line = rstrip(raw_line);
        if (line.empty() || line.substr(0, 2) == "##")
            continue;

        // Header line
        if (line.substr(0, 4) == "#CHR") {
            header = split(line, '\t');
            std::vector<std::string> inds_ID(header.begin() + 9, header.end());

            gt << "scaffold\tposition\teffect\tvartype\tflag\tref";
            for (const auto& id : inds_ID) gt << "\t" << id;
            gt << "\n";

            if (low_cov == "YES") {
                gt_cov << "scaffold\tposition\teffect\tvartype\tflag\tref";
                for (const auto& id : inds_ID) gt_cov << "\t" << id;
                gt_cov << "\n";
            }
            continue;
        }

        // Data line: replace phased '|' with unphased '/'
        line = replaceAll(line, "|", "/");
        std::vector<std::string> fields = split(line, '\t');

        // Extract genotype fields (FORMAT column 9 onward)
        std::vector<std::string> all_geno;
        for (size_t i = 9; i < fields.size(); ++i)
            all_geno.push_back(split(fields[i], ':')[0]);

        // Extract allelic depth for low-cov resampling (FORMAT subfield index 2)
        std::vector<std::string> all_coverage;
        if (low_cov == "YES") {
            for (size_t i = 9; i < fields.size(); ++i) {
                auto fmt = split(fields[i], ':');
                if (fmt.size() > 2) {
                    auto ad_parts = split(fmt[2], ',');
                    if (ad_parts.size() >= 2)
                        all_coverage.push_back(ad_parts[0] + ":" + ad_parts[1]);
                    else
                        all_coverage.push_back(".");
                } else {
                    all_coverage.push_back(".");
                }
            }
        }

        // Get alleles per population
        std::vector<std::string> alleles_pop1    = get_alleles_from_index(header, fields, pop1);
        std::vector<std::string> alleles_pop2    = get_alleles_from_index(header, fields, pop2);
        std::vector<std::string> alleles_pop_out = get_alleles_from_index(header, fields, pop_out);

        // Non-missing subsets (for missingness filter)
        std::vector<std::string> alleles_pop1_nomiss, alleles_pop2_nomiss, alleles_pop_out_nomiss;
        for (const auto& a : alleles_pop1)    if (a != ".") alleles_pop1_nomiss.push_back(a);
        for (const auto& a : alleles_pop2)    if (a != ".") alleles_pop2_nomiss.push_back(a);
        for (const auto& a : alleles_pop_out) if (a != ".") alleles_pop_out_nomiss.push_back(a);

        // Determine ancestral allele based on polarization strategy
        int anc_all = 0;
        std::string flag;

        std::vector<std::string> empty;
        if (polX == "POP_OUT") {
            auto res = get_ancestral_allele(alleles_pop_out, alleles_pop1, alleles_pop2);
            anc_all = res.first; flag = res.second;
        } else if (polX == "POP1") {
            auto res = get_ancestral_allele(empty, alleles_pop1, empty);
            anc_all = res.first; flag = res.second;
        } else if (polX == "POP2") {
            auto res = get_ancestral_allele(empty, empty, alleles_pop2);
            anc_all = res.first; flag = res.second;
        } else if (polX == "POP1POP2") {
            auto res = get_ancestral_allele(empty, alleles_pop1, alleles_pop2);
            anc_all = res.first; flag = res.second;
        } else { // REF or default: no repolarization
            anc_all = 0;
            flag    = "reference";
        }

        // Apply missingness filter
        bool pass_miss =
            (int)alleles_pop_out_nomiss.size() > (m0 * 2) - 1 &&
            (int)alleles_pop1_nomiss.size()    > (m1 * 2) - 1 &&
            (int)alleles_pop2_nomiss.size()    > (m2 * 2) - 1;

        if (!pass_miss) {
            ++missing_loci;
            continue;
        }

        // Extract locus metadata
        const std::string& scaffold    = fields[0];
        const std::string& position    = fields[1];
        const std::string& vcf_ref_all = fields[3];
        const std::string& vcf_alt_all = fields[4];
        const std::string& info_field  = fields[7];

        // Only process loci with SnpEff annotation
        if (!contains(info_field, "ANN"))
            continue;

        ++valid_loci;

        // Parse SnpEff ANN field
        std::string ann_section = info_field.substr(info_field.find("ANN=") + 4);
        auto ann_parts = split(ann_section, '/');
        std::string annotation = (ann_parts.size() > 1) ? ann_parts[1] : "";
        std::string effect     = (ann_parts.size() > 2) ? ann_parts[2] : "";

        std::string vartype;
        if      (contains(annotation, "missense"))   vartype = "missense";
        else if (contains(annotation, "synonymous")) vartype = "synonymous";
        else if (contains(annotation, "intergenic")) vartype = "intergenic";
        else if (contains(annotation, "intron"))     vartype = "intron";
        else if (contains(annotation, "downstream")) vartype = "downstream";
        else if (contains(annotation, "upstream"))   vartype = "upstream";
        else                                          vartype = "else";

        // Write locus prefix to output
        gt     << scaffold << "\t" << position << "\t" << effect << "\t"
               << vartype  << "\t" << flag;
        if (low_cov == "YES")
            gt_cov << scaffold << "\t" << position << "\t" << effect << "\t"
                   << vartype  << "\t" << flag;

        // Normalise missing genotypes ('.' -> './.')
        std::vector<std::string> clean_geno;
        for (const auto& g : all_geno)
            clean_geno.push_back((g == ".") ? "./." : g);

        // Repolarize if ancestral allele is the ALT allele
        std::vector<std::string> new_geno;
        if (anc_all == 0) {
            new_geno = clean_geno;
            gt << "\t" << vcf_ref_all << "\t";
            if (low_cov == "YES") gt_cov << "\t" << vcf_ref_all << "\t";
        } else {
            for (const auto& g : clean_geno) {
                std::string inv = replace_geno(g);
                inv = replaceAll(inv, "|", "/");
                new_geno.push_back(inv);
            }
            gt << "\t" << vcf_alt_all << "\t";
            if (low_cov == "YES") gt_cov << "\t" << vcf_alt_all << "\t";
            ++inverted_loci;
        }

        // Encode genotypes as derived-allele dosage and write
        for (size_t i = 0; i < new_geno.size(); ++i) {
            gt << replace_geno_coding(new_geno[i]);
            if (i + 1 < new_geno.size()) gt << "\t";
        }
        gt << "\n";

        // Low-coverage: write one randomly sampled allele per individual
        if (low_cov == "YES") {
            auto allele_random = random_allele_sampling(all_coverage, anc_all, rng);
            for (size_t i = 0; i < allele_random.size(); ++i) {
                gt_cov << allele_random[i];
                if (i + 1 < allele_random.size()) gt_cov << "\t";
            }
            gt_cov << "\n";
        }
    }

    gt.close();
    if (low_cov == "YES") gt_cov.close();

    std::cout << "\n" << valid_loci << " loci recorded to file, of which "
              << inverted_loci << " repolarized to ancestral, "
              << missing_loci  << " not passing missingness filter\n";
}

// =============================================================================
//  main — command-line argument parser
// =============================================================================
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <vcf_file> [options]\n"
                  << "Options:\n"
                  << "  --p1     <file>         Pop1 sample list\n"
                  << "  --p2     <file>         Pop2 sample list\n"
                  << "  --p0     <file>         Outgroup sample list\n"
                  << "  --m1     <int>          Min individuals Pop1 (default: 0)\n"
                  << "  --m2     <int>          Min individuals Pop2 (default: 0)\n"
                  << "  --m0     <int>          Min individuals outgroup (default: 0)\n"
                  << "  --polX   <string>       Polarization strategy: REF|POP_OUT|POP1|POP2|POP1POP2 (default: REF)\n"
                  << "  --low_cov <YES|NO>      Low-coverage resampling (default: NO)\n";
        return 1;
    }

    std::string vcf_file = argv[1];
    std::string p1, p2, p0;
    int m1 = 0, m2 = 0, m0 = 0;
    std::string polX    = "REF";
    std::string low_cov = "NO";

    for (int i = 2; i < argc; i += 2) {
        std::string key = argv[i];
        if (i + 1 >= argc) {
            std::cerr << "Missing value for argument: " << key << "\n";
            return 1;
        }
        std::string val = argv[i + 1];
        if      (key == "--p1")      p1      = val;
        else if (key == "--p2")      p2      = val;
        else if (key == "--p0")      p0      = val;
        else if (key == "--m1")      m1      = std::stoi(val);
        else if (key == "--m2")      m2      = std::stoi(val);
        else if (key == "--m0")      m0      = std::stoi(val);
        else if (key == "--polX")    polX    = val;
        else if (key == "--low_cov") low_cov = val;
        else {
            std::cerr << "Unknown argument: " << key << "\n";
            return 1;
        }
    }

    try {
        write_gt(vcf_file, p1, p2, p0, m1, m2, m0, polX, low_cov);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

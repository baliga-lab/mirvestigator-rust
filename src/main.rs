extern crate ftp;
extern crate flate2;
extern crate rayon;
extern crate clap;

use std::io;
use std::io::Read;
use std::collections::HashMap;
use std::fs::File;
use json;
use json::JsonValue;
use rayon::prelude::*;
//use clap::{Arg, App};
use clap::Parser;

const KMER8: &str = "8mer";
const KMER7_M8: &str = "7mer-m8";
const KMER7_A1: &str = "7mer-a1";
const NA: &str = "NA";

/*
An implementation of miRvestigator in Rust.
*/
fn all_kmers(alphabet: &[char], k: u32, seqs: Vec<String>, depth: u32) -> Vec<String> {
    if depth == 0 {
        let mut res = Vec::new();
        for c in alphabet {
            let mut s = String::new();
            s.push(*c);
            res.push(s);
        }
        all_kmers(alphabet, k, res, depth + 1)
    } else {
        if depth == k {
            seqs
        }
        else {
            let mut res = Vec::new();
            for s in seqs {
                for c in alphabet {
                    let mut snew = s.clone();
                    snew.push(*c);
                    res.push(snew);
                }
            }
            all_kmers(alphabet, k, res, depth + 1)
        }
    }
}

fn complement(seq: &str, comp_dict: &HashMap<char, char>) -> String {
    let mut s = String::new();
    for c in seq.chars() {
        match comp_dict.get(&c) {
            Some(comp) => s.push(*comp),
            None => panic!("no complement for {}", c)
        }
    }
    s
}

fn make_comp_dict() -> HashMap<char, char> {
    let mut comp_dict = HashMap::new();
    comp_dict.insert('A', 'T');
    comp_dict.insert('T', 'A');
    comp_dict.insert('G', 'C');
    comp_dict.insert('C', 'G');
    comp_dict.insert('U', 'A');
    comp_dict
}

fn is_valid_mirna(s: &str, minor: bool, p5: bool, p3: bool) -> bool {
    let mirna_id = String::from(s);
    (minor || !mirna_id.contains('*')) && (p5 || !mirna_id.contains("-5p"))
        && (p3 || !mirna_id.contains("-3p"))
}

fn read_mirnas(path: &str, seed_start: usize, seed_end: usize, minor: bool,
               p5: bool, p3: bool) -> HashMap<String, String> {
    use std::io::BufRead;
    use flate2::read::GzDecoder;

    println!("PATH: {}", path);
    let file = match File::open(path) {
        Err(why) => panic!("can't open mirbase input file: {}", why),
        Ok(file) => file,
    };
    let gz = GzDecoder::new(file);
    let mut reader = io::BufReader::new(gz);
    let mut mirna = String::new();
    let mut seq = String::new();
    let mut mi_rnas = HashMap::new();
    let comp_dict = make_comp_dict();

    loop {
        // description line
        let num_bytes = reader.read_line(&mut mirna).unwrap();
        if num_bytes == 0 {
            break;
        }
        // sequence line
        reader.read_line(&mut seq).unwrap();
        if mirna.starts_with(">hsa") {
            let mirna_comps: Vec<&str>  = mirna.trim_start_matches('>').split(' ').collect();
            let mirna_id = mirna_comps[0];
            if is_valid_mirna(&mirna_id, minor, p5, p3) {
                unsafe {
                    // reverse complement
                    let s = seq.trim().get_unchecked(seed_start..seed_end);
                    let rcs = complement(s.chars().rev().collect::<String>().as_str(),
                                         &comp_dict);
                    mi_rnas.insert(String::from(mirna_id), String::from(rcs));
                }
            }
        }
        // clear the buffers
        mirna.clear();
        seq.clear();
    }
    let mut mirna_uniq: HashMap<String, Vec<String>> = HashMap::new();
    for (mirna, seq) in mi_rnas {
        if !mirna_uniq.contains_key(&seq) {
            let mut v = Vec::new();
            v.push(mirna);
            mirna_uniq.insert(seq, v);
        } else {
            match mirna_uniq.get_mut(&seq) {
                Some(v) => v.push(mirna),
                None => panic!("should not happen")
            }
        }
    }
    let mut mirna_merged: HashMap<String, String> = HashMap::new();
    for (seed, elems) in mirna_uniq {
        let s = elems.join("_");
        mirna_merged.insert(s, seed);
    }
    mirna_merged
}

fn trim_seqs<'a>(mirnas: &'a HashMap<String, String>,
                 start: usize, stop: usize) -> HashMap<&'a String, String> {
    let mut res: HashMap<&String, String> = HashMap::new();
    for (mirna, seq) in mirnas {
        unsafe {
            let s = seq.get_unchecked(start..stop);
            res.insert(mirna, s.to_string());
        }
    }
    res
}

fn read_seqs(path: &str) -> String {
    use std::io::BufRead;

    let mut line = String::new();
    let mut result = String::new();
    let f = match File::open(path) {
        Err(why) => panic!("open error: {}", why),
        Ok(f) => f,
    };
    let mut reader = io::BufReader::new(f);
    let mut i = 0;
    loop {
        let num_bytes = reader.read_line(&mut line).unwrap();
        if num_bytes == 0 {
            break;
        }
        if i > 0 {
            result.push('_');
        }
        result.push_str(line.trim());
        line.clear();
        i += 1;
    }
    result
}

fn get_kmers(alphabet: &[char], k:u32, p3utr_seqs: &String, verbose: bool) -> Vec<String> {
    if verbose {
        eprintln!("Screening out {}mers not present in 3' UTRs...", k);
    }
    let km = all_kmers(alphabet, k, Vec::new(), 0);
    // filter by the p3utr sequences
    let mut kmers = Vec::new();
    for s in km {
        if p3utr_seqs.contains(s.as_str()) {
            kmers.push(s);
        }
    }
    kmers
}

fn add_hmm_state<'a>(states: &mut Vec<&'a String>, sp: &mut HashMap<&'a String, f64>,
                     state: &'a String, prob: f64) {
    states.push(state);
    sp.insert(state, prob);
}

fn add_hmm_transition<'a>(tp: &mut HashMap<&'a String, HashMap<&'a String, f64>>,
                      state1: &'a String, state2: &'a String, prob: f64) {
    if !tp.contains_key(state1) {
        tp.insert(state1, HashMap::new());
    }
    tp.get_mut(state1).unwrap().insert(state2, prob);
}

fn add_hmm_emission<'a>(ep: &mut HashMap<&'a String, HashMap<char, f64>>,
                        state1: &'a String, base: char, prob: f64) {
    if !ep.contains_key(state1) {
        ep.insert(state1, HashMap::new());
    }
    ep.get_mut(state1).unwrap().insert(base, prob);
}

fn print_viterbi_result(max_prob: f64, vit_mirnas: &Vec<&String>, label: &str) {
    eprint!("max viterbi probability ({}): {} [", label, max_prob);
    for m in vit_mirnas {
        eprint!(" {}", m);
    }
    eprintln!("]");
}

fn viterbi_kmers(kmers: &Vec<String>, states: &Vec<&String>,
                 sp: &HashMap<&String, f64>, tp: &HashMap<&String, HashMap<&String, f64>>,
                 ep: &HashMap<&String, HashMap<char, f64>>,
                 max_prob: f64, label: &str, verbose: bool) -> u32 {
    let mut hits = 0;
    let mut vitps = Vec::new();
    for kmer in kmers {
        let obs = kmer.chars().collect::<Vec<char>>();
        let prob = viterbi(&obs, states, sp, tp, ep);
        if prob > max_prob {
            hits = 2;
            break;
        } else if prob == max_prob {
            hits += 1;
            if hits > 1 {
                break;
            }
        }
        vitps.push(prob);
    }
    if hits <= 1 && verbose {
        eprintln!("=> {} match !", label);
    }
    hits
}

/* We borrow the results values from the key parameter's keys */
fn viterbi_mirnas<'a>(key_mirnas: &'a HashMap<&String, String>,
                      val_mirnas: &HashMap<&String, String>,
                      states: &Vec<&String>,
                      sp: &HashMap<&String, f64>, tp: &HashMap<&String, HashMap<&String, f64>>,
                      ep: &HashMap<&String, HashMap<char, f64>>) -> (f64, Vec<&'a String>) {
    let mut max_prob: f64 = -1.0;
    let mut vit_mirnas: Vec<&String> = Vec::new();
    for mirna in key_mirnas.keys() {
        let obs = val_mirnas.get(mirna).unwrap().chars().collect::<Vec<char>>();
        let prob = viterbi(&obs, &states, &sp, &tp, &ep);
        if prob > max_prob {
            max_prob = prob;
            vit_mirnas.clear();
            vit_mirnas.push(mirna);
        } else if prob == max_prob {
            vit_mirnas.push(mirna);
        }
    }
    (max_prob, vit_mirnas)
}

fn detect_pssm<'a, 'b>(pssm: &'a JsonValue,
               mirnas_6mer: &'b HashMap<&String, String>,
               mirnas_7mer_m8: &HashMap<&String, String>,
               mirnas_7mer_a1: &HashMap<&String, String>,
               mirnas_8mer: &HashMap<&String, String>,
               kmers7: &Vec<String>,
               kmers8: &Vec<String>,
               wobble: bool, wobble_cut: f64,
               verbose: bool) -> Vec<(&'a str, Vec<&'b String>, &'static str)> {
    eprintln!("Building HMM model for {}", pssm["name"]);
    let matrix = &pssm["matrix"];
    let max_pssm_i = matrix.len();
    let nm1 = String::from("NM1");
    let nm2 = String::from("NM2");
    // holder for indexed state names
    let mut wobble_states = Vec::new();
    let mut pssm_states = Vec::new();

    // list of HMM states
    let mut states = Vec::new();
    // state probabilities
    let mut sp = HashMap::new();
    // transition probabilities, which is a map from
    // State -> State -> probability
    let mut tp: HashMap<&String, HashMap<&String, f64>> = HashMap::new();

    // emission probabilities, which is a map from
    // State -> Base -> probability
    let mut ep: HashMap<&String, HashMap<char, f64>> = HashMap::new();

    // initialize state names
    for i in 0..max_pssm_i {
        pssm_states.push(format!("PSSM{}", i));
        wobble_states.push(format!("WOBBLE{}", i));
    }

    add_hmm_state(&mut states, &mut sp, &nm1, 1.0 / (max_pssm_i + 1) as f64);
    add_hmm_state(&mut states, &mut sp, &nm2, 0.0);

    // add PSSM states
    for i in 0..max_pssm_i {
        add_hmm_state(&mut states, &mut sp, &pssm_states[i],
                      1.0 / (max_pssm_i + 1) as f64);
        if wobble {
            add_hmm_state(&mut states, &mut sp, &wobble_states[i], 0.0);
        }
    }

    // add transitions
    // NM1
    let nm1_2_nm1 = 0.01;
    let left_over1 = (1.0 - nm1_2_nm1) / max_pssm_i as f64;

    add_hmm_transition(&mut tp, &nm1, &nm1, nm1_2_nm1);
    add_hmm_transition(&mut tp, &nm1, &nm2, 0.0);

    for i in 0..max_pssm_i {
        add_hmm_transition(&mut tp, &nm1, &pssm_states[i], left_over1);
        if wobble {
            add_hmm_transition(&mut tp, &nm1, &wobble_states[i], 0.0);
        }
    }
    // NM2
    add_hmm_transition(&mut tp, &nm2, &nm1, 0.0);
    add_hmm_transition(&mut tp, &nm2, &nm2, 1.0);
    for i in 0..max_pssm_i {
        add_hmm_transition(&mut tp, &nm2, &pssm_states[i], 0.0);
        if wobble {
            add_hmm_transition(&mut tp, &nm2, &wobble_states[i], 0.0);
        }
    }
    // PSSMis
    for i in 0..max_pssm_i {
        add_hmm_transition(&mut tp, &pssm_states[i], &nm1, 0.0);
        add_hmm_transition(&mut tp, &pssm_states[i], &nm2, 0.01);
        if wobble {
            add_hmm_transition(&mut tp, &wobble_states[i], &nm1, 0.0);
            add_hmm_transition(&mut tp, &wobble_states[i], &nm2, 0.01);
        }
        if i == (max_pssm_i - 1) {
            add_hmm_transition(&mut tp, &pssm_states[i], &nm2, 1.0);
            add_hmm_transition(&mut tp, &wobble_states[i], &nm2, 1.0);
        }

        for j in 0..max_pssm_i {
            if j == i + 1 {
                if wobble {
                    // allow wobbly matches if T is >= wobble_cut
                    if matrix[i+1][2].as_f64().unwrap() >= wobble_cut ||
                        matrix[i+1][3].as_f64().unwrap() >= wobble_cut {
                            add_hmm_transition(&mut tp, &pssm_states[i], &pssm_states[j], 0.8);
                            add_hmm_transition(&mut tp, &pssm_states[i], &wobble_states[j], 0.19);
                        } else {
                            // otherwise don't allow wobbly matches
                            add_hmm_transition(&mut tp, &pssm_states[i], &pssm_states[j], 0.99);
                            add_hmm_transition(&mut tp, &pssm_states[i], &wobble_states[j], 0.0);
                        }
                    add_hmm_transition(&mut tp, &wobble_states[i], &pssm_states[j], 1.0);
                    add_hmm_transition(&mut tp, &wobble_states[i], &wobble_states[j], 0.0);
                } else {
                    add_hmm_transition(&mut tp, &pssm_states[i], &pssm_states[j], 0.99);
                }
            } else {
                add_hmm_transition(&mut tp, &pssm_states[i], &pssm_states[j], 0.0);
                if wobble {
                    add_hmm_transition(&mut tp, &pssm_states[i], &wobble_states[j], 0.0);
                    add_hmm_transition(&mut tp, &wobble_states[i], &pssm_states[j], 0.0);
                    add_hmm_transition(&mut tp, &wobble_states[i], &wobble_states[j], 0.0);
                }
            }
        }
    }
    // Emission probabilities
    add_hmm_emission(&mut ep, &nm1, 'A', 0.25);
    add_hmm_emission(&mut ep, &nm1, 'C', 0.25);
    add_hmm_emission(&mut ep, &nm1, 'G', 0.25);
    add_hmm_emission(&mut ep, &nm1, 'T', 0.25);

    add_hmm_emission(&mut ep, &nm2, 'A', 0.25);
    add_hmm_emission(&mut ep, &nm2, 'C', 0.25);
    add_hmm_emission(&mut ep, &nm2, 'G', 0.25);
    add_hmm_emission(&mut ep, &nm2, 'T', 0.25);

    for i in 0..max_pssm_i {
        add_hmm_emission(&mut ep, &pssm_states[i], 'A', matrix[i][0].as_f64().unwrap());
        add_hmm_emission(&mut ep, &pssm_states[i], 'C', matrix[i][1].as_f64().unwrap());
        add_hmm_emission(&mut ep, &pssm_states[i], 'G', matrix[i][2].as_f64().unwrap());
        add_hmm_emission(&mut ep, &pssm_states[i], 'T', matrix[i][3].as_f64().unwrap());
        if wobble {
            // if motif has both G and T/U probability >= wobble_cut or random (0.25)
            if matrix[i][2].as_f64().unwrap() >= wobble_cut &&
                matrix[i][3].as_f64().unwrap() >= wobble_cut {
                    add_hmm_emission(&mut ep, &wobble_states[i], 'A', 0.5);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'C', 0.5);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'G', 0.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'T', 0.0);
                } else if matrix[i][2].as_f64().unwrap() >= wobble_cut {
                    // if motif has G >= wobble_cut or random (0.25)
                    add_hmm_emission(&mut ep, &wobble_states[i], 'A', 1.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'C', 0.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'G', 0.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'T', 0.0);
                } else if matrix[i][3].as_f64().unwrap() >= wobble_cut {
                    // if motif has U >= wobble_cut or random (0.25)
                    add_hmm_emission(&mut ep, &wobble_states[i], 'A', 0.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'C', 1.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'G', 0.0);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'T', 0.0);
                } else {
                    add_hmm_emission(&mut ep, &wobble_states[i], 'A', 0.25);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'C', 0.25);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'G', 0.25);
                    add_hmm_emission(&mut ep, &wobble_states[i], 'T', 0.25);
                }
        }
    }

    // Do Viterbi for all miRNAs
    let mut matches: Vec<(&str, Vec<&String>, &str)> = Vec::new();
    let (max_prob, vit_mirnas) = viterbi_mirnas(mirnas_6mer, mirnas_8mer,
                                                &states, &sp, &tp, &ep);
    if verbose {
        print_viterbi_result(max_prob, &vit_mirnas, "all");
    }
    let hits_8mer = viterbi_kmers(kmers8, &states, &sp, &tp, &ep, max_prob, KMER8, verbose);
    if hits_8mer <= 1 {
        // add to output (PSSMname,[mirnas],<label>)
        matches.push((pssm["name"].as_str().unwrap(), vit_mirnas, KMER8));
    }

    if hits_8mer > 1 {
        // try for perfect 7mer-m8
        let (max_prob, vit_mirnas) = viterbi_mirnas(mirnas_6mer, mirnas_7mer_m8,
                                                    &states, &sp, &tp, &ep);
        if verbose {
            print_viterbi_result(max_prob, &vit_mirnas, "7mer-m8");
        }
        let hits_7mer_m8 = viterbi_kmers(kmers7, &states, &sp, &tp, &ep, max_prob, KMER7_M8, verbose);
        if hits_7mer_m8 <= 1 {
            // add to output (PSSMname,[mirnas],<label>)
            matches.push((pssm["name"].as_str().unwrap(), vit_mirnas, KMER7_M8));
        }

        // finally try for perfect 7mer-a1
        let (max_prob, vit_mirnas) = viterbi_mirnas(mirnas_6mer, mirnas_7mer_a1,
                                                    &states, &sp, &tp, &ep);
        if verbose {
            print_viterbi_result(max_prob, &vit_mirnas, "7mer-a1");
        }
        let hits_7mer_a1 = viterbi_kmers(kmers7, &states, &sp, &tp, &ep, max_prob, KMER7_A1, verbose);
        if hits_7mer_a1 <= 1 {
            // add to output (PSSMname,[mirnas],<label>)
            matches.push((pssm["name"].as_str().unwrap(), vit_mirnas, KMER7_A1));
        }

        if hits_8mer > 1 && hits_7mer_m8 > 1 && hits_7mer_a1 > 1 {
            if verbose {
                eprintln!("=> no match !");
            }
            // Write to output (PSSM,[NA],NA)
            let empty = Vec::new();
            matches.push((pssm["name"].as_str().unwrap(), empty, NA));
        }
    }
    matches
}

/*
 * miRvestigator uses only the probability result of the Viterbi algorithm,
 * so we only run the actually used parts of it
 */
fn viterbi(obs: &Vec<char>, states: &Vec<&String>, start_p: &HashMap<&String, f64>,
           trans_p: &HashMap<&String, HashMap<&String, f64>>,
           emit_p: &HashMap<&String, HashMap<char, f64>>) -> f64 {
    // initialize v for T(0)
    let mut v = Vec::new();
    let mut ymap = HashMap::new();
    for y in states {
        // for now, we keep that ugly unwrap() cascade until we can verify this stuff works
        ymap.insert(y, start_p.get(y).unwrap() * emit_p.get(y).unwrap().get(&obs[0]).unwrap());
    }
    v.push(ymap);

    for t in 1..obs.len() {
        v.push(HashMap::new());
        // push a new time point
        for y in states {
            let mut max_prob = -1.0;
            for y0 in states {
                let tmp = v[t-1].get(y0).unwrap() * trans_p.get(y0).unwrap().get(y).unwrap()
                    * emit_p.get(y).unwrap().get(&obs[t]).unwrap();
                if tmp > max_prob {
                    max_prob = tmp;
                }
            }
            v[t].insert(y, max_prob);
        }
    }
    // return the max probability at time point t(n)
    let mut max_prob = -1.0;
    for y in states {
        let tmp: f64 = *v[obs.len()-1].get(y).unwrap();
        if tmp > max_prob {
            max_prob = tmp;
        }
    }
    max_prob
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about=None)]
struct Args {
    seqs: String,
    pssms: String,
    mirbase: String,
    outfile: String
}


fn main() -> io::Result<()> {
    use itertools::join;
    use std::io::Write;
    use std::sync::{Arc, Mutex};

    // Parameters
    let wobble = true;
    let wobble_cut = 0.25;
    let alphabet = ['A', 'C', 'G', 'T'];
    let minor = true;
    let p5 = true;
    let p3 = true;
    let verbose = true;

    let args = Args::parse();

    // p3utr_seqs are solely used for filtering the kmer lists
    let p3utr_seqs = read_seqs(&args.seqs);

    // Read in the JSON pssms file
    let mut pssms_json = String::new();
    let mut f2 = File::open(&args.pssms)?;
    let _ = f2.read_to_string(&mut pssms_json);
    let pssms = json::parse(pssms_json.as_str()).unwrap();
    if verbose {
        eprintln!("# seqs: {}", p3utr_seqs.len());
        eprintln!("# pssms: {}", pssms.len());
    }

    /* This HashMap is used to complement and reverse complement */
    let mirnas = read_mirnas(&args.mirbase, 0, 8, minor, p5, p3);
    let mirnas_6mer = trim_seqs(&mirnas, 0, 6);
    let mirnas_7mer_m8 = trim_seqs(&mirnas, 1, 8);
    let mirnas_7mer_a1 = trim_seqs(&mirnas, 0, 7);
    let mirnas_8mer = trim_seqs(&mirnas, 0, 8);

    // build a kmer list that is also present in p3utr_seqs
    let kmers7 = get_kmers(&alphabet, 7, &p3utr_seqs, verbose);
    let kmers8 = get_kmers(&alphabet, 8, &p3utr_seqs, verbose);

    if verbose {
        eprintln!("# final 7mers: {}, # final 8mers: {}", kmers7.len(), kmers8.len());
    }

    // build HMM model here
    let mut pssm_vec = Vec::new();
    for i in 0..pssms.len() {
        pssm_vec.push(&pssms[i]);
    }

    // TODO: collect all the matches together and write them out
    let output = Arc::new(Mutex::new(File::create(args.outfile)?));

    pssm_vec.par_iter().for_each(|pssm| {
        let matches = detect_pssm(pssm,
                                  &mirnas_6mer,
                                  &mirnas_7mer_m8, &mirnas_7mer_a1,
                                  &mirnas_8mer,
                                  &kmers7, &kmers8,
                                  wobble, wobble_cut, verbose);
        for (pssm_name, mirnas, label) in matches {
            // result to stdout
            let mirna_str = if mirnas.len() > 0 {
                join(mirnas, "_")
            } else {
                String::from(NA)
            };
            let _res = writeln!(&mut (output.lock().unwrap()), "{},{},{}", pssm_name, mirna_str, label);
        }
    });

    eprintln!("HMM Detection done.");
    Ok(())
}

// use std::env;
use std::fs::File;
use std::path::Path;
use std::env;
use std::io::prelude::*;

use std::io::{self, BufRead};


mod module_triangulation;
use crate::module_triangulation::TTriangulation;


fn parse_bool_token(s: &str) -> bool {
    let t = s.trim().to_ascii_lowercase();
    t == "true" || t == ".true." || t == "t"
}

fn ensure_size_2d<T: Clone>(mat: &mut Vec<Vec<T>>, rows: usize, cols: usize, fill: T) {
    while mat.len() < rows {
        mat.push(Vec::new());
    }

    for r in 0..rows {
        while mat[r].len() < cols {
            mat[r].push(fill.clone());
        }
    }
}

fn parse_ltri_entry(s: &str) -> io::Result<(usize, usize, i32)> {
    let open = s.find('(').ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Ltri sem '('"))?;
    let comma = s[open+1..].find(',').ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Ltri sem ','"))? + open + 1;
    let close = s[comma+1..].find(')').ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Ltri sem ')'"))? + comma + 1;
    let eq = s.find('=').ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Ltri sem '='"))?;

    let i_str = s[open+1..comma].trim();
    let j_str = s[comma+1..close].trim();
    let v_str = s[eq+1..].trim();

    let i = i_str.parse::<usize>().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "índice i inválido"))?;
    let j = j_str.parse::<usize>().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "índice j inválido"))?;
    let v = v_str.parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "valor v inválido"))?;

    Ok((i, j, v))
}


fn main() {

    let args: Vec<String> = env::args().collect();

    let case: Result<usize, _> = if args.len() > 1 {args[1].to_string().parse()} else {Ok(1)};
    let number_of_threads: Result<usize, _> = if args.len() > 2 {args[2].to_string().parse()} else {Ok(2)};
    let input_filename: Result<String, _> = if args.len() > 3 {args[3].to_string().parse()} else {Ok("1_000.txt".to_string())};
    
    // let intermediary_data_path = "../../../Data/Intermediary/";
    let intermediary_data_path = "./Data/Intermediary/";

    // let input_filename = "example_data.txt";
    // let input_filename = "input_data.txt";
    // let input_filename = "inicial.txt";
    // let input_filename = "inicial_1_000.txt";
    // let input_filename = "inicial_5_000.txt";
    // let input_filename = "inicial_10_000.txt";
    // let input_filename = "final_25_000.txt";
    // let input_filename = "final_50_000.txt";
    // let input_filename = "intermediario_75_000.txt";
    // let input_filename = "inicial_100_000.txt";

    let mut node_x: Vec<f32> = Vec::new();
    let mut node_y: Vec<f32> = Vec::new();
    let mut node_z: Vec<f32> = Vec::new();
    let mut number_of_nodes: usize = 0;

    let mut list: Vec<i32> = Vec::new();
    let mut lptr: Vec<i32> = Vec::new();
    let mut lend: Vec<i32> = Vec::new();
    let mut lnew: i32 = 0;
    let mut near: Vec<i32> = Vec::new();
    let mut next_tri: Vec<i32> = Vec::new();
    let mut ltri: Vec<Vec<i32>> = Vec::new();
    let mut bnodes: Vec<i32> = Vec::new();
    let mut dist: Vec<f32> = Vec::new();
    let mut number_of_triangles: usize = 0;
    let mut number_of_arcs: usize = 0;
    let mut number_of_boundary_nodes: usize = 0;
    let mut must_swap: bool = false;
    let mut swap_node: i32 = 0;
    let mut min_x: f32 = 0.;
    let mut max_x: f32 = 0.;
    let mut min_y: f32 = 0.;
    let mut max_y: f32 = 0.;
    
    let intermediary_data_path = Path::new(intermediary_data_path).join(&input_filename.clone().unwrap());
    let file = File::open(intermediary_data_path).expect("intermediary data file not found");
    let reader = io::BufReader::new(file);

    let mut section = String::new();
    
    for line in reader.lines() {
        let line = line.expect("Erro ao obter linha");
        let s = line.trim();
        if s.is_empty() {
            continue;
        }

        // Marcação de seções
        if s.starts_with("XT:")        { section = "XT".into();        continue; }
        if s.starts_with("YT:")        { section = "YT".into();        continue; }
        if s.starts_with("ZT:")        { section = "ZT".into();        continue; }
        if s.starts_with("List:")      { section = "List".into();      continue; }
        if s.starts_with("Lptr:")      { section = "Lptr".into();      continue; }
        if s.starts_with("Lend:")      { section = "Lend".into();      continue; }
        if s.starts_with("Near:")      { section = "Near".into();      continue; }
        if s.starts_with("NextTri:")   { section = "NextTri".into();   continue; }
        if s.starts_with("Ltri:")      { section = "Ltri".into();      continue; }
        if s.starts_with("BNodes:")    { section = "BNodes".into();    continue; }
        if s.starts_with("Dist:")      { section = "Dist".into();      continue; }

        // Escalares (podem vir combinados na mesma linha)
        if s.contains("number_of_nodes:") {
            if let Some(x) = s.split("number_of_nodes:").nth(1) {
                number_of_nodes = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("number_of_triangles:") {
            if let Some(x) = s.split("number_of_triangles:").nth(1) {
                number_of_triangles = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("number_of_arcs:") {
            if let Some(x) = s.split("number_of_arcs:").nth(1) {
                number_of_arcs = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("number_of_boundary_nodes:") {
            if let Some(x) = s.split("number_of_boundary_nodes:").nth(1) {
                number_of_boundary_nodes = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("Lnew:") {
            if let Some(x) = s.split("Lnew:").nth(1) {
                lnew = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("must_swap:") {
            if let Some(x) = s.split("must_swap:").nth(1) {
                let tok = x.trim().split_whitespace().next().unwrap_or("");
                must_swap = parse_bool_token(tok);
            }
        }
        if s.contains("swap_node:") {
            if let Some(x) = s.split("swap_node:").nth(1) {
                swap_node = x.trim().split_whitespace().next().unwrap().parse().unwrap();
            }
        }
        if s.contains("min_x:") && s.contains("max_x:") {
            if let Some(x) = s.split("min_x:").nth(1) {
                let mut it = x.trim().split_whitespace();
                min_x = it.next().unwrap().parse().unwrap();
            }
            if let Some(x) = s.split("max_x:").nth(1) {
                let mut it = x.trim().split_whitespace();
                max_x = it.next().unwrap().parse().unwrap();
            }
        }
        if s.contains("min_y:") && s.contains("max_y:") {
            if let Some(x) = s.split("min_y:").nth(1) {
                let mut it = x.trim().split_whitespace();
                min_y = it.next().unwrap().parse().unwrap();
            }
            if let Some(x) = s.split("max_y:").nth(1) {
                let mut it = x.trim().split_whitespace();
                max_y = it.next().unwrap().parse().unwrap();
            }
        }

        // Vetores / matriz
        match section.as_str() {
            "XT" | "YT" | "ZT" | "List" | "Lptr" | "Lend" | "Near" | "NextTri" | "BNodes" | "Dist" => {
                // Formato esperado: "Name(i) = valor"
                if let Some(rhs) = s.split('=').nth(1) {
                    let val_str = rhs.trim();
                    match section.as_str() {
                        "XT"      => node_x.push(val_str.parse().unwrap()),
                        "YT"      => node_y.push(val_str.parse().unwrap()),
                        "ZT"      => node_z.push(val_str.parse().unwrap()),
                        "List"    => list.push(val_str.parse().unwrap()),
                        "Lptr"    => lptr.push(val_str.parse().unwrap()),
                        "Lend"    => lend.push(val_str.parse().unwrap()),
                        "Near"    => near.push(val_str.parse().unwrap()),
                        "NextTri" => next_tri.push(val_str.parse().unwrap()),
                        "BNodes"  => bnodes.push(val_str.parse().unwrap()),
                        "Dist"    => dist.push(val_str.parse().unwrap()),
                        _ => {}
                    }
                }
            }
            "Ltri" => {
                let (i_idx, j_idx, value) = parse_ltri_entry(s).expect("Erro ao extrair informações");
                let r = i_idx - 1;
                let c = j_idx - 1;

                let needed_rows = r + 1;
                let needed_cols = c + 1;
                ensure_size_2d(&mut ltri, needed_rows, needed_cols, 0);

                ltri[r][c] = value;
            }
            _ => {}
        }
    }

    let _ = TTriangulation::construct_triangulation_xyz(
        case.unwrap(), // triangulation_id
        number_of_threads.unwrap(),
        number_of_nodes, // num_nodes
        node_x,
        node_y,
        node_z,
        0.0, // tolerance
        list,
        lptr,
        lend,
        lnew,
        near,
        next_tri,
        ltri,
        bnodes,
        dist,
        number_of_triangles,
        number_of_arcs,
        number_of_boundary_nodes,
        must_swap,
        swap_node,
        min_x,
        max_x,
        min_y,
        max_y,
        &mut 0, // stat
        &input_filename.unwrap()
    );
}
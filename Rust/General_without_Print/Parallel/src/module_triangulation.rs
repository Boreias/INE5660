use std::i32;
use std::time::Instant;
use std::fs::File;
use std::io::{Write, Result};
use std::path::Path;

use std::collections::HashMap;
use std::sync::Arc;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

#[derive(Debug, Clone)]
struct TDirectedEdge {
    edge_id: i32,
    start_node: i32,
    end_node: i32,
    voronoi_x: f32,
    voronoi_y: f32,
    length: f32,
    slope: f32,
    boundary: bool,
    counter_clock_edge: Option<Box<TDirectedEdge>>,
}

#[derive(Debug, Clone)]
struct TNode {
    x: f32,
    y: f32,
    z: f32,
    n_neighbor: i32,
    boundary: bool,
    state: i32,
    voronoi_area: f32,
    first_edge: Option<Box<TDirectedEdge>>,
    flow_exit: Option<Box<TDirectedEdge>>,
}

#[derive(Debug, Clone, Copy)]
struct TTriangle {
    center_x: f32,
    center_y: f32,
    radius: f32,
    area: f32,
    aspect_ratio: f32,
}

#[derive(Debug, Clone)]
struct TReach {
    edge: Option<Box<TDirectedEdge>>,
    n_strahler: i32,
    drainage_area: f32,
    next: Option<Box<TReach>>,
}

#[derive(Debug, Clone)]
pub struct TTriangulation {
    instance_id: usize,
    xt: Vec<f32>,
    yt: Vec<f32>,
    zt: Vec<f32>,
    list: Vec<i32>,
    lptr: Vec<i32>,
    lend: Vec<i32>,
    lnew: i32,
    near: Vec<i32>,
    next_tri: Vec<i32>,
    ltri: Vec<Vec<i32>>,
    bnodes: Vec<i32>,
    dist: Vec<f32>,
    number_of_nodes: usize,
    number_of_triangles: usize,
    number_of_arcs: usize,
    number_of_boundary_nodes: usize,
    must_swap: bool,
    swap_node: i32,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    nodes: Vec<TNode>,
    directed_edges: Vec<TDirectedEdge>,
    triangles: Vec<TTriangle>,
    have_height_values: bool,
    first_reach: Option<Box<TReach>>,
    n_reaches: i32,
}


impl TTriangulation {

    pub fn construct_triangulation_xyz(
        triangulation_id: usize,
        number_of_threads: usize,
        mut number_of_nodes: usize,
        mut node_x: Vec<f32>,
        mut node_y: Vec<f32>,
        node_z: Vec<f32>,
        tolerance: f32,
        list: Vec<i32>,
        lptr: Vec<i32>,
        lend: Vec<i32>,
        lnew: i32,
        near: Vec<i32>,
        next_tri: Vec<i32>,
        ltri: Vec<Vec<i32>>,
        bnodes: Vec<i32>,
        dist: Vec<f32>,
        number_of_triangles: usize,
        number_of_arcs: usize,
        number_of_boundary_nodes: usize,
        must_swap: bool,
        swap_node: i32,
        min_x: f32,
        max_x: f32,
        min_y: f32,
        max_y: f32,
        stat: &mut i32,
        filename: &str
    ) -> Result<()> {

        if number_of_nodes < 3 {
            println!("Insuficient number of triangulation nodes");
            return Ok(());
        }

        let mut triangulation = TTriangulation {
            instance_id: triangulation_id,
            xt: node_x.clone(),
            yt: node_y.clone(),
            zt: node_z.clone(),
            list: list,
            lptr: lptr,
            lend: lend,
            lnew: lnew,
            near: near,
            next_tri: next_tri,
            ltri: ltri,
            bnodes: bnodes,
            dist: dist,
            number_of_nodes,
            number_of_triangles: number_of_triangles,
            number_of_arcs: number_of_arcs,
            number_of_boundary_nodes: number_of_boundary_nodes,
            must_swap: must_swap,
            swap_node: swap_node,
            min_x: min_x,
            max_x: max_x,
            min_y: min_y,
            max_y: max_y,
            nodes: Vec::new(),
            directed_edges: Vec::new(),
            triangles: Vec::new(),
            have_height_values: false,
            first_reach: None,
            n_reaches: 0,
        };


        triangulation.construct_nodes(&node_x, &node_y, &node_z, number_of_threads);

        triangulation.construct_triangles(number_of_threads);

        triangulation.construct_directed_edges(number_of_threads);


        Ok(())
    }

    fn construct_directed_edges(&mut self, num_threads: usize) {
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("failed to build rayon pool");

        let number_of_triangles = self.number_of_triangles;
        let expected_candidates = (number_of_triangles as usize) * 6;

        #[derive(Clone)]
        struct EdgeCand {
            start: i32,
            end: i32,
            edge_id: i32,
            voronoi_x: f32,
            voronoi_y: f32,
            length: f32,
            slope: f32,
            boundary: bool,
        }

        let candidates: Vec<EdgeCand> = pool.install(|| {
            (0..number_of_triangles)
                .into_par_iter()
                .flat_map(|i_t| {
                    let v1 = self.ltri[0][i_t];
                    let v2 = self.ltri[1][i_t];
                    let v3 = self.ltri[2][i_t];

                    let n1 = self.ltri[3][i_t];
                    let n2 = self.ltri[4][i_t];
                    let n3 = self.ltri[5][i_t];

                    let e1 = self.ltri[6][i_t];
                    let e2 = self.ltri[7][i_t];
                    let e3 = self.ltri[8][i_t];

                    let tri_cx = self.triangles[i_t].center_x;
                    let tri_cy = self.triangles[i_t].center_y;

                    let calc_props = |start: i32, end: i32, nb: i32| -> (f32, f32, bool) {
                        let s = &self.nodes[(start - 1) as usize];
                        let e = &self.nodes[(end - 1) as usize];
                        let dx = s.x - e.x;
                        let dy = s.y - e.y;
                        let dz = s.z - e.z;
                        let length = (dx * dx + dy * dy).sqrt();
                        let slope = if length != 0.0 { dz / length } else { 0.0 };
                        let boundary_edge = s.boundary && e.boundary && nb == 0;
                        (length, slope, boundary_edge)
                    };

                    let mut local: Vec<EdgeCand> = Vec::with_capacity(6);

                    let candidates_for_triangle = [
                        (v2, v3, e1, n1),
                        (v3, v1, e2, n2),
                        (v1, v2, e3, n3),
                    ];

                    for &(start, end, edge_id, nb) in &candidates_for_triangle {
                        let (length, slope, boundary_edge) = calc_props(start, end, nb);

                        let (vx1, vy1) = if nb != 0 {
                            let nb_idx = (nb - 1) as usize;
                            (self.triangles[nb_idx].center_x, self.triangles[nb_idx].center_y)
                        } else {
                            (f32::NAN, f32::NAN)
                        };

                        local.push(EdgeCand {
                            start,
                            end,
                            edge_id,
                            voronoi_x: vx1,
                            voronoi_y: vy1,
                            length,
                            slope,
                            boundary: boundary_edge,
                        });

                        local.push(EdgeCand {
                            start: end,
                            end: start,
                            edge_id,
                            voronoi_x: tri_cx,
                            voronoi_y: tri_cy,
                            length,
                            slope: -slope,
                            boundary: boundary_edge,
                        });
                    }

                    local.into_par_iter()
                })
                .collect::<Vec<EdgeCand>>()
        });

        let mut edge_map: HashMap<(i32, i32), usize> = HashMap::with_capacity(candidates.len() / 2 + 1);
        let mut edges: Vec<TDirectedEdge> = Vec::with_capacity(candidates.len() / 2 + 1);

        for cand in candidates.into_iter() {
            let key = (cand.start, cand.end);
            if !edge_map.contains_key(&key) {
                let idx = edges.len();
                edges.push(TDirectedEdge {
                    edge_id: cand.edge_id,
                    start_node: cand.start,
                    end_node: cand.end,
                    voronoi_x: cand.voronoi_x,
                    voronoi_y: cand.voronoi_y,
                    length: cand.length,
                    slope: cand.slope,
                    boundary: cand.boundary,
                    counter_clock_edge: None,
                });
                edge_map.insert(key, idx);
            }
        }

        for idx in 0..edges.len() {
            if edges[idx].counter_clock_edge.is_some() {
                continue;
            }
            let twin_key = (edges[idx].end_node, edges[idx].start_node);
            if let Some(&twin_idx) = edge_map.get(&twin_key) {
                edges[idx].counter_clock_edge = Some(Box::new(edges[twin_idx].clone()));
                edges[twin_idx].counter_clock_edge = Some(Box::new(edges[idx].clone()));
            }
        }

        self.directed_edges = edges;

        use std::collections::HashMap as Map;
        let mut first_edge_for_node: Map<i32, usize> = Map::with_capacity(self.number_of_nodes);

        for (idx, de) in self.directed_edges.iter().enumerate() {
            first_edge_for_node.entry(de.start_node).or_insert(idx);
        }

        for i_n in 0..self.number_of_nodes {
            let node_id = (i_n as i32) + 1;
            if let Some(&idx) = first_edge_for_node.get(&node_id) {
                self.nodes[i_n].first_edge = Some(Box::new(self.directed_edges[idx].clone()));
            } else {
                self.nodes[i_n].first_edge = None;
            }
        }
    }

    fn construct_nodes(
        &mut self,
        node_x: &Vec<f32>,
        node_y: &Vec<f32>,
        node_z: &Vec<f32>,
        num_threads: usize,
    ) {
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        pool.install(|| {
            let results: Vec<(TNode, f32, f32, f32, f32)> = (0..self.number_of_nodes)
                .into_par_iter()
                .map(|i_n| {
                    let x = node_x[i_n];
                    let y = node_y[i_n];
                    let z = node_z[i_n];

                    let n_neighbor = self.nbcnt(self.lend[i_n] - 1);

                    let mut boundary = false;
                    for &b in &self.bnodes {
                        if b == (i_n as i32) {
                            boundary = true;
                            break;
                        }
                    }

                    (
                        TNode {
                            x,
                            y,
                            z,
                            n_neighbor,
                            boundary,
                            state: 0,
                            voronoi_area: 0.0,
                            first_edge: None,
                            flow_exit: None,
                        },
                        x,
                        x,
                        y,
                        y,
                    )
                })
                .collect();

            let mut nodes = Vec::with_capacity(self.number_of_nodes);
            let mut min_x = f32::MAX;
            let mut max_x = f32::MIN;
            let mut min_y = f32::MAX;
            let mut max_y = f32::MIN;

            for (node, minx, maxx, miny, maxy) in results {
                nodes.push(node);
                if minx < min_x { min_x = minx; }
                if maxx > max_x { max_x = maxx; }
                if miny < min_y { min_y = miny; }
                if maxy > max_y { max_y = maxy; }
            }

            self.nodes = nodes;
            self.min_x = min_x;
            self.max_x = max_x;
            self.min_y = min_y;
            self.max_y = max_y;
            self.have_height_values = true;
        });
    }

    fn construct_triangles(&mut self, num_threads: usize) {
        let pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("failed to build rayon thread pool");

        let triangles: Vec<TTriangle> = pool.install(|| {
            (0..self.number_of_triangles)
                .into_par_iter()
                .map(|i_t| {
                    let idx1 = (self.ltri[0][i_t] - 1) as usize;
                    let idx2 = (self.ltri[1][i_t] - 1) as usize;
                    let idx3 = (self.ltri[2][i_t] - 1) as usize;

                    let x1 = self.nodes[idx1].x;
                    let y1 = self.nodes[idx1].y;
                    let x2 = self.nodes[idx2].x;
                    let y2 = self.nodes[idx2].y;
                    let x3 = self.nodes[idx3].x;
                    let y3 = self.nodes[idx3].y;

                    let u1 = x3 - x2;
                    let u2 = x1 - x3;
                    let u3 = x2 - x1;

                    let v1 = y3 - y2;
                    let v2 = y1 - y3;
                    let v3 = y2 - y1;

                    let area = (u1 * v2 - u2 * v1) * 0.5_f32;

                    let (center_x, center_y, radius, aspect_ratio) =
                        if area.abs() > f32::EPSILON {
                            let ds1 = x1 * x1 + y1 * y1;
                            let ds2 = x2 * x2 + y2 * y2;
                            let ds3 = x3 * x3 + y3 * y3;

                            let fx = -(ds1 * v1 + ds2 * v2 + ds3 * v3);
                            let fy =  (ds1 * u1 + ds2 * u2 + ds3 * u3);

                            let cx = fx / (4.0_f32 * area);
                            let cy = fy / (4.0_f32 * area);

                            let r = ((cx - x1).powi(2) + (cy - y1).powi(2)).sqrt();

                            let e1 = u1 * u1 + v1 * v1;
                            let e2 = u2 * u2 + v2 * v2;
                            let e3 = u3 * u3 + v3 * v3;

                            let perimeter = e1.sqrt() + e2.sqrt() + e3.sqrt();

                            let ar = if r > 0.0 && perimeter > 0.0 {
                                2.0_f32 * area.abs() / (perimeter * r)
                            } else {
                                0.0_f32
                            };

                            (cx, cy, r, ar)
                        } else {
                            (0.0_f32, 0.0_f32, 0.0_f32, 0.0_f32)
                        };

                    TTriangle {
                        center_x,
                        center_y,
                        radius,
                        area,
                        aspect_ratio,
                    }
                })
                .collect()
        });

        self.triangles = triangles;
    }

    fn nbcnt (&self, lpl: i32) -> i32 {
        let mut lp = lpl;
        let mut k = 1;

        loop {
            lp = self.lptr[lp as usize] - 1;
            if lp == lpl {
                break;
            }
            k += 1;
        }

        k as i32
    }
}
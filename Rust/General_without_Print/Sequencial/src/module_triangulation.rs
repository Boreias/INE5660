use std::io::Result;

#[derive(Debug, Clone)]
struct TDirectedEdge {
    edge_id: usize,
    start_node: usize,
    end_node: usize,
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
    n_neighbor: usize,
    boundary: bool,
    state: usize,
    voronoi_area: f32,
    first_edge: Option<Box<TDirectedEdge>>,
    flow_exit: Option<Box<TDirectedEdge>>,
}

#[derive(Debug, Clone)]
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
    n_strahler: usize,
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
    lptr: Vec<usize>,
    lend: Vec<usize>,
    lnew: usize,
    near: Vec<usize>,
    next_tri: Vec<usize>,
    ltri: Vec<Vec<usize>>,
    bnodes: Vec<usize>,
    dist: Vec<f32>,
    number_of_nodes: usize,
    number_of_triangles: usize,
    number_of_arcs: usize,
    number_of_boundary_nodes: usize,
    must_swap: bool,
    swap_node: usize,
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    nodes: Vec<TNode>,
    directed_edges: Vec<TDirectedEdge>,
    triangles: Vec<TTriangle>,
    have_height_values: bool,
    first_reach: Option<Box<TReach>>,
    n_reaches: usize,
}


impl TTriangulation {

    pub fn construct_triangulation_xyz(
        triangulation_id: usize,
        number_of_nodes: usize,
        node_x: Vec<f32>,
        node_y: Vec<f32>,
        node_z: Vec<f32>,
        tolerance: f32,
        list: Vec<i32>,
        lptr: Vec<usize>,
        lend: Vec<usize>,
        lnew: usize,
        near: Vec<usize>,
        next_tri: Vec<usize>,
        ltri: Vec<Vec<usize>>,
        bnodes: Vec<usize>,
        dist: Vec<f32>,
        number_of_triangles: usize,
        number_of_arcs: usize,
        number_of_boundary_nodes: usize,
        must_swap: bool,
        swap_node: usize,
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


        triangulation.construct_nodes(&node_x, &node_y, &node_z);

        triangulation.construct_triangles();

        triangulation.construct_directed_edges();


        Ok(())
    }

    fn construct_directed_edges(&mut self) {

        self.directed_edges = vec![TDirectedEdge {
            edge_id: 0,
            start_node: 0,
            end_node: 0,
            voronoi_x: 0.0,
            voronoi_y: 0.0,
            length: 0.,
            slope: 0.,
            boundary: false,
            counter_clock_edge: None,
        }; 2 * self.number_of_arcs];

        let mut edge_so_far: usize = 0;

        for it in 0..self.number_of_triangles {
            let v1 = self.ltri[0][it];
            let v2 = self.ltri[1][it];
            let v3 = self.ltri[2][it];

            let n1 = self.ltri[3][it];
            let n2 = self.ltri[4][it];
            let n3 = self.ltri[5][it];

            let e1 = self.ltri[6][it];
            let e2 = self.ltri[7][it];
            let e3 = self.ltri[8][it];

            if !self.already_in_list(&edge_so_far, v2, v3) {
                self.add_edge(&mut edge_so_far, e1, v2, v3, n1, it);
            }

            if !self.already_in_list(&edge_so_far, v3, v1 ) {
                self.add_edge(&mut edge_so_far, e2, v3, v1, n2, it);
            }

            if !self.already_in_list(&edge_so_far, v1, v2) {
                self.add_edge(&mut edge_so_far, e3, v1, v2, n3, it);
            }
        }

        for ie in 0..self.number_of_arcs {
            let start_node = self.directed_edges[ie].start_node;
            let end_node = self.directed_edges[ie].end_node;

            let lpl = self.lend[(start_node - 1) as usize];
            let mut lp = lpl;
            let mut k = 0;

            let mut nd: i32 = 0;
            let mut nabor = vec![0; 1000];
            let mut end_node_counter_clock: i32 = 0;

            loop {
                lp = self.lptr[(lp - 1) as usize];
                nd = self.list[(lp - 1) as usize];
                nabor[k] = nd;
                if lp == lpl {
                    break;
                }
                k += 1;
            }
            k += 1;

            if nd <= 0 {
                nabor[k] = -nd;
            }

            for i in 1..=k {
                if nabor[i] == end_node as i32 {
                    if i == k {
                        end_node_counter_clock = nabor[1];
                    } else {
                        end_node_counter_clock = nabor[i + 1];
                    }
                    break;
                }
            }

            for ie2 in 0..self.number_of_arcs {
                if self.directed_edges[ie2].start_node == start_node {
                    if self.directed_edges[ie2].end_node as i32 == end_node_counter_clock {
                        self.directed_edges[ie].counter_clock_edge = Some(Box::new(self.directed_edges[ie2].clone()));
                        break;
                    }
                }
            }
        }

        for i_n in 0..self.number_of_nodes {
            for i_e in 0..self.number_of_arcs {
                if self.directed_edges[i_e].start_node == i_n {
                    self.nodes[i_n].first_edge = Some(Box::new(self.directed_edges[i_e].clone()));
                    break;
                }
            }
        }
    }

    fn already_in_list(&self, edge_so_far: &usize, start_node: usize, end_node: usize) -> bool {
        for i in (0..*edge_so_far).rev() {
            if self.directed_edges[i].start_node == start_node {
                if self.directed_edges[i].end_node == end_node {
                    return true;
                }
            }
        }

        return false;
    }

    fn add_edge(&mut self, edge_so_far: &mut usize, edge_id: usize, start_node: usize, end_node: usize, nb: usize, it: usize) {

        let dx = self.nodes[start_node - 1].x - self.nodes[end_node - 1].x;
        let dy = self.nodes[start_node - 1].y - self.nodes[end_node - 1].y;
        let dz = self.nodes[start_node - 1].z - self.nodes[end_node - 1].z;

        let length = (dx.powf(2.0) + dy.powf(2.0)).sqrt();
        let slope = dz / length;

        let mut boundary_edge = false;

        if self.nodes[start_node - 1].boundary && self.nodes[end_node - 1].boundary && nb == 0 {
            boundary_edge = true;
        }

        self.directed_edges[*edge_so_far].edge_id = edge_id;
        self.directed_edges[*edge_so_far].start_node = start_node;
        self.directed_edges[*edge_so_far].end_node = end_node;
        self.directed_edges[*edge_so_far].length = length;
        self.directed_edges[*edge_so_far].slope = slope;
        self.directed_edges[*edge_so_far].boundary = boundary_edge;

        if nb != 0 {
            self.directed_edges[*edge_so_far].voronoi_x = self.triangles[nb - 1].center_x;
            self.directed_edges[*edge_so_far].voronoi_y = self.triangles[nb - 1].center_y;
        } else {
            self.directed_edges[*edge_so_far].voronoi_x = f32::NAN;
            self.directed_edges[*edge_so_far].voronoi_y = f32::NAN;
        }

        *edge_so_far += 1;

        self.directed_edges[*edge_so_far].edge_id = edge_id;
        self.directed_edges[*edge_so_far].start_node = end_node;
        self.directed_edges[*edge_so_far].end_node = start_node;
        self.directed_edges[*edge_so_far].voronoi_x = self.triangles[it].center_x;
        self.directed_edges[*edge_so_far].voronoi_y = self.triangles[it].center_y;
        self.directed_edges[*edge_so_far].length = length;
        self.directed_edges[*edge_so_far].slope = -slope;
        self.directed_edges[*edge_so_far].boundary = boundary_edge;

        *edge_so_far += 1;
    }

    fn construct_nodes (&mut self, node_x: &Vec<f32>, node_y: &Vec<f32>, node_z: &Vec<f32>) {
        self.nodes = vec![TNode {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            n_neighbor: 0,
            boundary: false,
            state: 0,
            voronoi_area: 0.0,
            first_edge: None,
            flow_exit: None,
        }; self.number_of_nodes];

        for i_n in 0..self.number_of_nodes {
            self.nodes[i_n].x = node_x[i_n];
            self.nodes[i_n].y = node_y[i_n];
            self.nodes[i_n].z = node_z[i_n];
            self.nodes[i_n].n_neighbor = self.nbcnt(self.lend[i_n] - 1);
            self.nodes[i_n].boundary = false;

            for i_b in 0..self.number_of_boundary_nodes {
                if self.bnodes[i_b] == i_n {
                    self.nodes[i_n].boundary = true;
                    break;
                }
            }

            if self.nodes[i_n].x < self.min_x { self.min_x = self.nodes[i_n].x; }
            if self.nodes[i_n].x > self.max_x { self.max_x = self.nodes[i_n].x; }
            if self.nodes[i_n].y < self.min_y { self.min_y = self.nodes[i_n].y; }
            if self.nodes[i_n].y > self.max_y { self.max_y = self.nodes[i_n].y; }
        }

        self.have_height_values = true;
    }

    fn construct_triangles (&mut self) {
        self.triangles = vec![TTriangle {
            center_x: 0.0,
            center_y: 0.0,
            radius: 0.0,
            area: 0.0,
            aspect_ratio: 0.0,
        }; self.number_of_triangles];

        for i_t in 0..self.number_of_triangles {
            self.circum(
                self.nodes[(self.ltri[0][i_t] - 1) as usize].x,
                self.nodes[(self.ltri[0][i_t] - 1) as usize].y,
                self.nodes[(self.ltri[1][i_t] - 1) as usize].x,
                self.nodes[(self.ltri[1][i_t] - 1) as usize].y,
                self.nodes[(self.ltri[2][i_t] - 1) as usize].x,
                self.nodes[(self.ltri[2][i_t] - 1) as usize].y,
                true,
                i_t
            );
        }
    }

    fn circum (
        &mut self,
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        x3: f32,
        y3: f32,
        ratio: bool,
        i_t: usize
    ) {
        let u: Vec<f32> = vec![x3 - x2, x1 - x3, x2 - x1];
        let v: Vec<f32> = vec![y3 - y2, y1 - y3, y2 - y1];

        self.triangles[i_t].area = (u[0] * v[1] - u[1] * v[0]) / 2.0;

        if self.triangles[i_t].area == 0.0 {
            if ratio {
                self.triangles[i_t].aspect_ratio = 0.0;
            }
            return;
        }

        let ds: Vec<f32> = vec![
            x1 * x1 + y1 * y1,
            x2 * x2 + y2 * y2,
            x3 * x3 + y3 * y3,
        ];

        let fx = -ds.iter().zip(&v).map(|(d, v)| d * v).sum::<f32>();
        let fy = ds.iter().zip(&u).map(|(d, u)| d * u).sum::<f32>();

        self.triangles[i_t].center_x = fx / (4.0 * self.triangles[i_t].area);
        self.triangles[i_t].center_y = fy / (4.0 * self.triangles[i_t].area);
        self.triangles[i_t].radius = ((self.triangles[i_t].center_x - x1).powi(2) + (self.triangles[i_t].center_y - y1).powi(2)).sqrt();

        if !ratio {
            return;
        }

        let ds: Vec<f32> = u.iter().zip(v.iter()).map(|(&ui, &vi)| ui.powi(2) + vi.powi(2)).collect();

        self.triangles[i_t].aspect_ratio = 2.0 * self.triangles[i_t].area.abs() / ((ds[0].sqrt() + ds[1].sqrt() + ds[2].sqrt()) * self.triangles[i_t].radius);
    }

    fn nbcnt (&self, lpl: usize) -> usize {
        let mut lp = lpl;
        let mut k = 1;

        loop {
            lp = self.lptr[lp] - 1;
            if lp == lpl {
                break;
            }
            k += 1;
        }

        k
    }
}
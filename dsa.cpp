#include "dsa.h"

void Heap::Swap(int k, int l){
    PII tmp = h[k]; h[k] = h[l]; h[l] = tmp;
	table[h[k].first] = k; table[h[l].first] = l;
}

void Heap::Up(int k) {
    while(k > 1){
        int pa = k / 2;
        if (h[k].second < h[pa].second){
            Swap(k, pa);
            k = pa;
        }
        else break;
    }
}

void Heap::Down(int k) {
    while(k * 2 <= tail){
        int son = k * 2;
        if(son < tail && h[son+1].second < h[son].second) son++;
        if (h[k].second > h[son].second){
            Swap(k, son);
            k = son;
        }
        else break;
    }
}

void Heap::Insert(PII val) {
    table[val.first] = ++tail;
    h[tail] = val;
    Up(tail);
}

void Heap::Modify(int k, PII val) {
    h[k] = val;
	Up(k);
	Down(k);
}

int Heap::Pop() {
    int top = h[1].first; st[top] = true;
    Swap(1, tail); tail--;
    Down(1);
    return top;
}

void DualGraph::Assign_Weight(Graph& graph) {
    double tot_Ang_Dist = 0.0, tot_Geod = 0.0;
    int pair_num = 0;

    // traverse all edges
    for(int i = 0; i < graph.edges.size(); ++i) {
        for(int j = i + 1; j < graph.edges.size(); ++j) {
            // find neighbour faces
            if ((graph.edges[i].v_id[0] == graph.edges[j].v_id[0] && graph.edges[i].v_id[1] == graph.edges[j].v_id[1]) || 
                (graph.edges[i].v_id[0] == graph.edges[j].v_id[1] && graph.edges[i].v_id[1] == graph.edges[j].v_id[0])){
                
                double Angle = Normal_Angle(graph.faces[graph.edges[i].f_id], graph.faces[graph.edges[j].f_id]);
                double Ang_Dist = Cal_Ang_Dist(graph.faces[graph.edges[i].f_id], graph.faces[graph.edges[j].f_id]);
                double Geod = Cal_Geo_Dist(graph.faces[graph.edges[i].f_id], graph.faces[graph.edges[j].f_id], graph.vertices[graph.edges[i].v_id[0]], graph.vertices[graph.edges[i].v_id[1]]);
                
                tot_Ang_Dist += Ang_Dist;
                tot_Geod += Geod;
                pair_num += 1;

				graph.faces[graph.edges[i].f_id].neighbours.push_back(Neighbour(graph.edges[j].f_id, Ang_Dist, Geod, Angle));
				graph.faces[graph.edges[j].f_id].neighbours.push_back(Neighbour(graph.edges[i].f_id, Ang_Dist, Geod, Angle));
            }
        }
    }

    graph.avg_Ang_Dist = tot_Ang_Dist / pair_num;
	graph.avg_Geod = tot_Geod / pair_num;
    // calculate weights of dual graph
	for(int i = 0; i < graph.faces.size(); ++i) {
        for(int j = 0; j < graph.faces[i].neighbours.size(); ++j) {
            graph.faces[i].neighbours[j].dist = DELTA * graph.faces[i].neighbours[j].Geod / graph.avg_Geod + 
                (1 - DELTA) * graph.faces[i].neighbours[j].Ang_Dist / graph.avg_Ang_Dist;
        }
    }
}

void DualGraph::Shortest_Path(Graph& graph) {
    // Dijkstra + Binary Heap
    int n = graph.faces.size(); // dual vertices (graph faces)
    w.resize(n); dis.resize(n); ne.resize(n); heap.Resize(n);

    // initialization
    for(int i = 0; i < n; ++i) {
        dis[i].resize(n);
        for(int j = 0; j < n; ++j) {
            dis[i][j] = j==i ? 0.0 : INF;
        }
        for(int j = 0; j < graph.faces[i].neighbours.size(); ++j) {
            w[i].push_back(graph.faces[i].neighbours[j].dist);
            ne[i].push_back(graph.faces[i].neighbours[j].f_id);
        }
    }

    // dijkstra, source: u
    for(int u = 0; u < n; ++u) {
        for (int i = 0; i < n; ++i) {
            heap.table[i] = 0;
            heap.st[i] = false;
        }
        heap.tail = 0;
        heap.Insert(std::make_pair(u, 0.0));
        while(heap.tail > 0){
            int top = heap.Pop();
            for(int i = 0; i < ne[top].size(); ++i){
                int v = ne[top][i];
                double d = w[top][i];
                if (!heap.st[v] && dis[u][v] > dis[u][top] + d) {
                    // relax
                    dis[u][v] = dis[u][top] + d;
                    if(heap.table[v] == 0) heap.Insert(std::make_pair(v, dis[u][v]));
                    else heap.Modify(heap.table[v], std::make_pair(v, dis[u][v]));
                }
            }
        }
    }
    std::cout << "shortest path algorithm finish!" << std::endl;
}

void GraphSolver::Init(Graph* _g, DualGraph* _dg, int l, std::vector<int>& s) {
    g = _g; dg = _dg; level = l;
    int n = g->faces.size();
    P_A.resize(n); P_B.resize(n); flow_net.resize(n+2);
    S = n; T = n + 1;
    if(l == 0){
        for (int i = 0; i < n; ++i) sub_v.push_back(i);
        sub_n = sub_v.size();
    }
    else{
        sub_v = s; sub_n = sub_v.size();
    }
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < g->faces[i].neighbours.size(); ++j) {
            FlowEdge fe;
            fe.id = g->faces[i].neighbours[j].f_id;
			fe.cap = 1.0 / (1.0 + g->faces[i].neighbours[j].Ang_Dist / g->avg_Ang_Dist);
			fe.flow = 0.0;
			flow_net[i].push_back(fe);
        }
    }
    avg_dist = Cal_Avg_Dist();
    if(level == 0) GLOBAL_AVG_DIST = avg_dist;
    Find_Max_Dist(); // init RepA, RepB
}

double GraphSolver::Cal_Avg_Dist() {
    double sum_dist = 0.0;
    for(int i = 0; i < sub_n; ++i) {
        for(int j = i + 1; j < sub_n; ++j) {
            sum_dist += dg->dis[sub_v[i]][sub_v[j]];
        }
    }
    double avg_dist = 2 * sum_dist / (sub_n * (sub_n - 1));
    return avg_dist;
}

void GraphSolver::Find_Max_Dist() {
    double max_dist = -INF;
    for(int i = 0; i < sub_n; ++i) {
        for(int j = i + 1; j < sub_n; ++j) {
            if(max_dist < dg->dis[sub_v[i]][sub_v[j]]) {
                REP_A = sub_v[i]; REP_B = sub_v[j];
                max_dist = dg->dis[sub_v[i]][sub_v[j]];
            }
        }
    }
    if(level == 0) GLOBAL_MAX_DIST = max_dist;
}

void GraphSolver::Solve() {
    int iter = 0;
    while(iter < MAX_ITER){
        int old_REPA = REP_A, old_REPB = REP_B;
        Assign_Prob();
        Update_Rep();
        if(REP_A == old_REPA && REP_B == old_REPB) break;
        iter++;
    }
    Assign_Label();
    Graph_Cut();
    LABEL_BASE += K; // ensure assigning different labels
    DFS_Solve();
}

void GraphSolver::Assign_Prob() {
    for(int i = 0; i < sub_n; ++i) {
        if(sub_v[i] == REP_A) {P_A[i] = 1.0; P_B[i] = 0.0;}
        else if(sub_v[i] == REP_B) {P_B[i] = 1.0; P_A[i] = 0.0;}
        else{
            P_A[i] = dg->dis[sub_v[i]][REP_B] / (dg->dis[sub_v[i]][REP_A] + dg->dis[sub_v[i]][REP_B]);
            P_B[i] = 1.0 - P_A[i];
        }
    }
}

void GraphSolver::Assign_Label() {
    int label, fuzzy_cnt = 0;
    int cnt[K] = {0};
    for(int i = 0; i < sub_n; ++i) {
        if(P_A[i] > P_B[i]) {label = A_FLAG;}
        else {label = B_FLAG;}
        if(abs(P_A[i] - P_B[i]) > PROB_THR) {cnt[label] += 1;}
        else {label = FUZZY; fuzzy_cnt += 1;}
        g->faces[sub_v[i]].label = LABEL_BASE + label;
    }
    std::cout << "level: " << level << " labelA: " << cnt[0] << ", labelB: " << cnt[1] << ", fuzzy: " << fuzzy_cnt << std::endl;
}

void GraphSolver::Update_Rep() {
    double F_A = INF, F_B = INF;
    double sum;
    for(int i = 0; i < sub_n; ++i){
        sum = 0.0;
        for(int j = 0; j < sub_n; ++j){
            sum += P_A[j] * dg->dis[sub_v[i]][sub_v[j]];
        }
        if(sum < F_A){
            F_A = sum; REP_A = sub_v[i];
        }

        sum = 0.0;
        for(int j = 0; j < sub_n; ++j){
            sum += P_B[j] * dg->dis[sub_v[i]][sub_v[j]];
        }
        if(sum < F_B){
            F_B = sum; REP_B = sub_v[i];
        }
    }
}

void GraphSolver::Init_Flow(std::vector<int>& p) {
    int n = g->faces.size();
    flow_net[S].clear();
    for(int i = 0; i < n; ++i) {
        if(p[i] == A_FLAG + 1){ // Set V_CA
            FlowEdge fe(i, INF, 0.0);
            flow_net[S].push_back(fe);
        }
        else if(p[i] == B_FLAG + 1){ // Set V_CB
            FlowEdge fe(T, INF, 0.0);
            flow_net[i].push_back(fe);
        }
    }
    p.push_back(OTHER_FLAG); // S
    p.push_back(OTHER_FLAG); // T
}

void GraphSolver::Max_Flow(std::vector<int>& p) {
    // Edmonds-Karp
    int n = g->faces.size();
    std::vector<FlowState> q;
    std::vector<bool> vis(n+2); // consider S, T
    double flow_sum = 0.0;
    int head = 0, tail = 0;
    for(;;) {
        head = tail = 0; q.clear();
        for(int i = 0; i < vis.size(); ++i) vis[i] = false;
        FlowState fs(S, 0, INF); // add S state
        q.push_back(fs); vis[S] = true;
        while(head <= tail && !vis[T]) { // BFS
            FlowState top = q[head++];
            for(int i = 0; i < flow_net[top.id].size(); ++i){
                int ne = flow_net[top.id][i].id;
                double inc = flow_net[top.id][i].cap - flow_net[top.id][i].flow;
                if (p[ne] == 0) continue; // not in the sub_graph
                if (vis[ne] || inc < EPS) continue;
                FlowState fs(ne, head - 1, std::min(top.inc, inc));
                q.push_back(fs); vis[ne] = true; tail++;
                if(ne == T) break; // stop until we find T
            }
        }
        if(!vis[T]) { // no more augmenting path
            for(int i = 1; i <= tail; ++i) p[q[i].id] = A_FLAG + 1;
            for(int i = 0; i < n; ++i) {if (p[i] == FUZZY) p[i] = B_FLAG + 1;} // eliminate fuzzy region
            break; // algorithm finish
        }
        double inc = q[tail].inc; flow_sum += inc;
        while (tail != 0){ // update flow
            int prev = q[tail].front;
            int src = q[prev].id;
			int dst = q[tail].id;
            for(int i = 0; i < flow_net[src].size(); ++i){
                if(flow_net[src][i].id == dst) flow_net[src][i].flow += inc; // forward flow
            }
            for(int i = 0; i < flow_net[dst].size(); ++i){
                if(flow_net[dst][i].id == src) flow_net[dst][i].flow -= inc; // backward flow
            }
			tail = prev;
		}
    }
    std::cout << "level: " << level << " max flow computed: " << flow_sum << std::endl;
}

void GraphSolver::Graph_Cut() {
    int n = g->faces.size();
    std::vector<int> partition; partition.resize(n);
    for(int i = 0; i < n; ++i) partition[i] = 0;
    for(int i = 0; i < sub_n; ++i){
        if(g->faces[sub_v[i]].label == LABEL_BASE + FUZZY){
            partition[sub_v[i]] = FUZZY; // Set V_C
            for(int j = 0; j < g->faces[sub_v[i]].neighbours.size(); ++j){
                int ne = g->faces[sub_v[i]].neighbours[j].f_id;
                if(g->faces[ne].label == LABEL_BASE + A_FLAG) partition[ne] = A_FLAG + 1; // Set V_CA
                else if(g->faces[ne].label == LABEL_BASE + B_FLAG) partition[ne] = B_FLAG + 1; // Set V_CB
            }
        }
    }
    Init_Flow(partition);
    Max_Flow(partition);
    for(int i = 0; i < sub_n; ++i){ // assign labels
        if(partition[sub_v[i]] == A_FLAG + 1) g->faces[sub_v[i]].label = LABEL_BASE + A_FLAG; // A
        else if(partition[sub_v[i]] == B_FLAG + 1) g->faces[sub_v[i]].label = LABEL_BASE + B_FLAG; // B
    }
}

void GraphSolver::DFS_Solve() {
    double rep_dist = dg->dis[REP_A][REP_B];
	if(level == MAX_DEPTH || (rep_dist / GLOBAL_MAX_DIST < REP_DIST_RATIO)) return;
    GraphSolver* gs_A = new GraphSolver();
    GraphSolver* gs_B = new GraphSolver();
    std::vector<int> sub_gA, sub_gB;
    for(int i = 0; i < sub_n; ++i) {
        if ((g->faces[sub_v[i]].label) % K == A_FLAG) sub_gA.push_back(sub_v[i]); // A
        if ((g->faces[sub_v[i]].label) % K == B_FLAG) sub_gB.push_back(sub_v[i]); // B
    }
    gs_A->Init(g, dg, level + 1, sub_gA);
    gs_B->Init(g, dg, level + 1, sub_gB);
	if(gs_A->avg_dist / GLOBAL_AVG_DIST > AVG_DIST_RATIO) gs_A->Solve();
    if(gs_B->avg_dist / GLOBAL_AVG_DIST > AVG_DIST_RATIO) gs_B->Solve();
}
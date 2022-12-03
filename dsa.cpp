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

void GraphSolver::Init(Graph* _g, DualGraph* _dg, int l) {
    g = _g; dg = _dg; level = l;
    int n = g->faces.size();
    P_A.resize(n); P_B.resize(n); flow_net.resize(n+2);
    S = n; T = n + 1;
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < g->faces[i].neighbours.size(); ++j) {
            FlowEdge fe;
            fe.id = g->faces[i].neighbours[j].f_id;
			fe.cap = 1.0 / (1.0 + g->faces[i].neighbours[j].Ang_Dist / g->avg_Ang_Dist);
			fe.flow = 0.0;
			flow_net[i].push_back(fe);
        }
    }
    Find_Max_Dist(); // init RepA, RepB
}

double GraphSolver::Find_Max_Dist() {
    double max_dist = -INF;
    int n = g->faces.size();
    for(int i = 0; i < n; ++i) {
        for(int j = i + 1; j < n; ++j) {
            if(max_dist < dg->dis[i][j]) {
                REP_A = i; REP_B = j;
                max_dist = dg->dis[i][j];
            }
        }
    }
    return max_dist;
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
    std::cout << "mesh decomposition solved!" << std::endl;
}

void GraphSolver::Assign_Prob() {
    int n = g->faces.size();
    for(int i = 0; i < n; ++i) {
        if(i == REP_A) {P_A[i] = 1.0; P_B[i] = 0.0;}
        else if(i == REP_B) {P_B[i] = 1.0; P_A[i] = 0.0;}
        else{
            P_A[i] = dg->dis[i][REP_B] / (dg->dis[i][REP_A] + dg->dis[i][REP_B]);
            P_B[i] = 1.0 - P_A[i];
        }
    }
}

void GraphSolver::Assign_Label() {
    int n = g->faces.size();
    int label, fuzzy_cnt = 0;
    int cnt[2] = {0, 0};
    for(int i = 0; i < n; ++i) {
        if(P_A[i] > P_B[i]) {label = 0;}
        else {label = 1;}
        if(abs(P_A[i] - P_B[i]) > PROB_THR) {cnt[label] += 1;}
        else {label = FUZZY; fuzzy_cnt += 1;}
        g->faces[i].label = label;
    }
    std::cout << "labelA: " << cnt[0] << ", labelB: " << cnt[1] << ", fuzzy: " << fuzzy_cnt << std::endl;
}

void GraphSolver::Update_Rep() {
    double F_A = INF, F_B = INF;
    double sum;
    int n = g->faces.size();
    for(int i = 0; i < n; ++i){
        sum = 0.0;
        for(int j = 0; j < n; ++j){
            sum += P_A[j] * dg->dis[i][j];
        }
        if(sum < F_A){
            F_A = sum; REP_A = i;
        }

        sum = 0.0;
        for(int j = 0; j < n; ++j){
            sum += P_B[j] * dg->dis[i][j];
        }
        if(sum < F_B){
            F_B = sum; REP_B = i;
        }
    }
}

void GraphSolver::Init_Flow(std::vector<int>& p) {
    int n = g->faces.size();
    flow_net[S].clear();
    for(int i = 0; i < n; ++i) {
        if(p[i] == 1){ // Set V_CA
            FlowEdge fe(i, INF, 0.0);
            flow_net[S].push_back(fe);
        }
        else if(p[i] == 2){ // Set V_CB
            if(flow_net[i][flow_net[i].size()-1].id != T){
                FlowEdge fe(T, INF, 0.0);
                flow_net[i].push_back(fe);
            }
        }
        else{ // Set C & others
			if (flow_net[i][flow_net[i].size()-1].id == T)
				flow_net[i].pop_back();
        }
    }
    p.push_back(4); // S
    p.push_back(4); // T
    for(int i = 0; i < n + 2; ++i){
        if(p[i] != 0) {
            for(int j = 0; j < flow_net[i].size(); ++j) flow_net[i][j].flow = 0.0;
        }
    }
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
                if (p[ne] == 0) continue;
                if (vis[ne] || inc < EPS) continue;
                FlowState fs(ne, head - 1, std::min(top.inc, inc));
                q.push_back(fs); vis[ne] = true; tail++;
                if(ne == T) break;
            }
        }
        if(!vis[T]) { // no more augmenting path
            for(int i = 1; i <= tail; ++i) p[q[i].id] = 1;
            for(int i = 0; i < n; ++i) {if (p[i] == 3) p[i] = 2;} // eliminate fuzzy region
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
    std::cout << "max flow computed: " << flow_sum << std::endl;
}

void GraphSolver::Graph_Cut() {
    int n = g->faces.size();
    std::vector<int> partition; partition.resize(n);
    for(int i = 0; i < n; ++i) partition[i] = 0;
    for(int i = 0; i < n; ++i){
        if(g->faces[i].label == FUZZY){
            partition[i] = 3; // Set V_C
            for(int j = 0; j < g->faces[i].neighbours.size(); ++j){
                int ne = g->faces[i].neighbours[j].f_id;
                if(g->faces[ne].label == 0) partition[ne] = 1; // Set V_CA
                else if(g->faces[ne].label == 1) partition[ne] = 2; // Set V_CB
            }
        }
    }
    Init_Flow(partition);
    Max_Flow(partition);
    for(int i = 0; i < n; ++i){ // assign labels
        if(partition[i] == 1) g->faces[i].label = 0; // A
        else if(partition[i] == 2) g->faces[i].label = 1; // B
    }
}
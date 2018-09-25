import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

class Solver6 {
	 
	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;
 
	private static final double START_TEMP = 100;
	private static final double END_TEMP = 10;
 
	private final Xor128 rand = new Xor128();
	private final long start = System.currentTimeMillis();
	private final long timeLimit;
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;
	private double temp;
 
	public Solver6(Data data) {
 
		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;
		this.timeLimit = data.timeLimit;
 
		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}
 
	public Result solve() {
		State bestState = new State();
		;
		loop: for (int len = 2; len >= 1; len--) {
			State state = new State();
			int L = (int) (Math.sqrt(Vemb) + 0.000001);
			Set<Integer> set = new HashSet<>();
			for (int i = 1; i <= V; i++) {
				set.add(i);
			}
 
			int[][] p = new int[V][2];
			for (int i = 0; i < V; i++) {
				p[i][0] = i + 1;
				p[i][1] = g[i + 1].length;
			}
			int LL = (int) (Math.sqrt(V) + 1.0000001);
			Arrays.sort(p, (o1, o2) -> o2[1] - o1[1]);
 
			int x = 0;
			int y = 0;
			for (int i = 0; i < V; i++) {
				while (state.map[x + y * L + 1] > 0) {
					x++;
					if (x >= LL) {
						y += 1;
						x = 0;
					}
					if (x + y * L + 1 > Vemb) {
						continue loop;
					}
				}
 
				boolean diag = false;
				if (x % len == 1 && y % len == 1 && x < L - len + 1 && y < L - len + 1) {
					diag = true;
				}
 
				List<Integer> list = new ArrayList<>();
				int now = x + y * L + 1;
				list.add(now);
				if (diag) {
					int q = now + L + 1;
					for (int l = 0; l < len - 1; l++) {
						if (q > Vemb) {
							continue loop;
						}
						list.add(q);
						q += L + 1;
					}
				} else {
					diag = false;
				}
 
				int max = -1;
				int maxV = 0;
				for (int j = 0; j < V; j++) {
					int v = diag ? p[j][0] : p[V - j - 1][0];
					if (!set.contains(v)) {
						continue;
					}
					int z = 0;
					for (int vEmb : list) {
						for (int uEmb : h[vEmb]) {
							int u = state.map[uEmb];
							z += d[u][v] ? 1 : 0;
						}
					}
					if (z > max) {
						maxV = v;
						max = z;
					}
				}
				for (int vEmb : list) {
					if (!state.set(maxV, vEmb)) {
						System.out.println("NG");
					}
				}
				set.remove(maxV);
			}
 
			if (bestState.score() < state.score()) {
				bestState = state;
			}
			break;
		}
 
		bestState = annealing(bestState);
 
		Result result = getResult(bestState, bestState.nodes.size, bestState.edgeCount);
		return result;
	}
 
	private State annealing(State state) {
		State bestState = new State(state);
		int bestScore = bestState.score();
		int currentScore = state.score();
		int L = (int) (Math.sqrt(Vemb) + 0.000001);
		while (!isEnd()) {
			int uEmb = state.nodes.array[rand.next() % state.nodes.size];
			int vEmb = state.nodes.array[rand.next() % state.nodes.size];
			;
			if (!swap(state, uEmb, vEmb, L)) {
				continue;
			}
			int nextScore = state.score();
			int scoreDiff = nextScore - currentScore;
 
			double probability = Math.exp(((double) scoreDiff) / temp);
			boolean forceNext = probability > (double) (rand.next() % Xor128.RAND_INF) / Xor128.RAND_INF;
 
			if (scoreDiff > 0 || forceNext) {
				currentScore = nextScore;
				if (bestScore < currentScore) {
					bestState = new State(state);
					bestScore = currentScore;
				}
			} else {
				swap(state, uEmb, vEmb, L);
			}
		}
 
		return bestState;
	}
 
	private boolean swap(State state, int uEmb, int vEmb, int L) {
		int v = state.map[vEmb];
		int u = state.map[uEmb];
		if (v == 0 || u == 0 || u == v) {
			return false;
		}
 
		int ulen = 0;
		int vlen = 0;
		for (; vEmb > L + 1 && state.map[vEmb - L - 1] == v; vEmb -= L + 1) {
		}
		for (; uEmb > L + 1 && state.map[uEmb - L - 1] == u; uEmb -= L + 1) {
		}
		for (int i = vEmb; i <= Vemb; i += L + 1) {
			if (state.map[i] == v) {
				state.remove(i);
				vlen++;
			} else {
				break;
			}
		}
		for (int i = uEmb; i <= Vemb; i += L + 1) {
			if (state.map[i] == u) {
				state.remove(i);
				ulen++;
			} else {
				break;
			}
		}
 
		for (int i = 0; i < ulen; i++) {
			if (!state.set(v, uEmb + i * (L + 1))) {
				System.out.println("NG");
			}
		}
		for (int i = 0; i < vlen; i++) {
			if (!state.set(u, vEmb + i * (L + 1))) {
				System.out.println("NG");
			}
		}
		return true;
	}
 
	private Result getResult(State state, int... params) {
		Result result = new Result(state.map, state.score(), time, loop, params);
		return result;
	}
 
	private boolean isEnd() {
		loop++;
		if (loop % 100 != 0) {
			return false;
		}
 
		long end = System.currentTimeMillis();
		time = end - start;
		temp = START_TEMP + (END_TEMP - START_TEMP) * time / timeLimit;
		return time >= timeLimit;
	}
 
	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;
 
		// vEmb->v
		public final int[] map;
 
		private final IntSet nodes;
		private final int[][] dc;
 
		private int edgeCount = 0;
 
		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			nodes = new IntSet(Vemb + 1);
		}
 
		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}
 
		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				set(v, vEmb);
			}
		}
 
		public boolean set(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			map[vEmb] = v;
			nodes.add(vEmb);
 
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}
 
		public int remove(int vEmb) {
			if (map[vEmb] == 0) {
				return 0;
			}
			int v = map[vEmb];
			map[vEmb] = 0;
			nodes.remove(vEmb);
 
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				dc[u][v]--;
				dc[v][u]--;
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount--;
				}
			}
			return v;
		}
 
		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
 
		public boolean isFull() {
			return edgeCount == E;
		}
	}
}

class Solver5 {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;


	private final Xor128 rand = new Xor128();
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;

	public Solver5(Data data) {

		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;

		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}

	public Result solve() {
		State bestState = new State();;
		loop: for (int len = 1; len < 5; len ++) {
			State state = new State();
			int L = (int)(Math.sqrt(Vemb) + 0.000001);
			Set<Integer> set = new HashSet<>();
			for (int i = 1; i <= V; i ++) {
				set.add(i);
			}
			
			int[][] p = new int[V][2];
			for (int i = 0; i < V; i ++) {
				p[i][0] = i + 1;
				p[i][1] = g[i + 1].length;
			}
			int LL = (int)(Math.sqrt(V) + 1.0000001);
			Arrays.sort(p, (o1, o2) -> o2[1] - o1[1]);

			int x = 0;
			int y = 0;
			for (int i = 0; i < V; i ++) {
				while (state.map[x + y * L + 1] > 0) {
					x ++;
					if (x >= LL) {
						y += 1;
						x = 0;
					}
					if (x + y * L + 1 > Vemb) {
						continue loop;
					}
				}
	
				boolean diag = false;
				if (x % len == 1 && y % len == 1 && x < L - len + 1 && y < L - len + 1) {
					diag = true;
				}
				
				List<Integer> list = new ArrayList<>();
				int now = x + y * L + 1;
				list.add(now);
				if (diag) {
					int q = now + L + 1;
					for (int l = 0; l < len - 1; l ++) {
						if (q > Vemb) {
							continue loop;
						}
						list.add(q);
						q += L + 1;
					}
				} else {
					diag = false;
				}
				
				int max = -1;
				int maxV = 0;
				for (int j = 0; j < V; j ++) {
					int v = diag ? p[j][0] : p[V - j - 1][0];
					if (!set.contains(v)) {
						continue;
					}
					int z = 0;
					for (int vEmb : list) {
						for (int uEmb : h[vEmb]) {
							int u = state.map[uEmb];
							z += d[u][v] ? 1 : 0;
						}
					}
					if (z > max) {
						maxV = v;
						max = z;
					}
				}
				for (int vEmb : list) {
					if (!state.addRoot(maxV, vEmb)) {
						System.out.println("NG");
					}
				}
				set.remove(maxV);
			}
			
			if (bestState.score() < state.score()) {
				bestState = state;
			}
		}

		Result result = getResult(bestState, bestState.nodes.size, bestState.edgeCount);
		return result;
	}

	private Result getResult(State state, int... params) {
		Result result = new Result(state.map, state.score(), time, loop, params);
		return result;
	}



	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;

		// vEmb->v
		public final int[] map;

		private final IntSet nodes;
		private final int[][] dc;

		private int edgeCount = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			nodes = new IntSet(Vemb + 1);
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}

		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				addRoot(v, vEmb);
			}
		}

		public boolean addRoot(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			map[vEmb] = v;
			nodes.add(vEmb);

			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}


		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
		
		public boolean isFull() {
			return edgeCount == E;
		}
	}
}



class Solver4 {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;

	private static final int ADD_RETRY = 10;

	private final Xor128 rand = new Xor128();
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;

	public Solver4(Data data) {

		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;

		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}

	public Result solve() {
		State state = new State();
		int L = (int)(Math.sqrt(Vemb) + 0.000001);
		Set<Integer> set = new HashSet<>();
		for (int i = 1; i <= V; i ++) {
			set.add(i);
		}
		boolean flg = Vemb > V * 1.5;
		
		int[][] p = new int[V][2];
		for (int i = 0; i < V; i ++) {
			p[i][0] = i + 1;
			p[i][1] = g[i + 1].length;
		}
		Arrays.sort(p, (o1, o2) -> o2[1] - o1[1]);
		
		int now = 0;
		for (int i = 0; i < V; i ++) {
			while (state.map[now + 1] > 0) {
				now ++;
			}

			int x = now % L;
			int y = now / L;
			boolean diag = false;
			if (x % 2 == 1 && y % 2 == 1 && x < L - 1) {
				diag = true;
			}
			
			List<Integer> list = new ArrayList<>();
			list.add(now + 1);
			if (flg && diag) {
				list.add(now + L + 2);
			} else {
				diag = false;
			}
			
			int max = -1;
			int maxV = 0;
			for (int j = 0; j < V; j ++) {
				int v = diag ? p[j][0] : p[V - j - 1][0];
				if (!set.contains(v)) {
					continue;
				}
				int z = 0;
				for (int vEmb : list) {
					for (int uEmb : h[vEmb]) {
						int u = state.map[uEmb];
						z += d[u][v] ? 1 : 0;
					}
				}
				if (z > max) {
					maxV = v;
					max = z;
				}
			}
			for (int vEmb : list) {
				state.addRoot(maxV, vEmb);
			}
			set.remove(maxV);
		}

		Result result = getResult(state, state.nodes.size, state.edgeCount);
		return result;
	}

	private Result getResult(State state, int... params) {
		Result result = new Result(state.map, state.score(), time, loop, params);
		return result;
	}


	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;

		// vEmb->v
		public final int[] map;

		private final IntSet leafs;
		private final IntSet nodes;
		private final int[] embbedEdgeCount;
		private final int[][] dc;

		private int edgeCount = 0;
		private int lastEmb = 0;
		private int lastV = 0;
		private int lastOperation = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			embbedEdgeCount = new int[Vemb + 1];
			leafs = new IntSet(Vemb + 1);
			nodes = new IntSet(Vemb + 1);
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			embbedEdgeCount = Arrays.copyOf(state.embbedEdgeCount, Vemb + 1);
			leafs = new IntSet(state.leafs);
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}

		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				addRoot(v, vEmb);
			}
		}

		public void addRoot(int v, int vEmb) {
			map[vEmb] = v;
			nodes.add(vEmb);

			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			
		}

		public boolean randomNext() {
			for (int k = 0; k < ADD_RETRY; k ++) {
				boolean ret = false;
				if (rand.next() % Xor128.RAND_INF < Xor128.RAND_INF) {
					ret = addRandomLeaf();
				} else {
					ret = removeRandomLeaf();
				}
				if (ret) {
					return true;
				}
			}
			return false;
		}

		public void undo() {
			if (lastOperation == OP_ADD) {
				removeLeaf(lastEmb);
			} else if (lastOperation == OP_REMOVE) {
				addLeaf(lastV, lastEmb);
			}
			lastOperation = 0;
		}

		public boolean addRandomLeaf() {
			if (nodes.size == 0) {
				return false;
			}
			int vEmb = nodes.array[rand.next() % nodes.size];
			int[] cand = h[vEmb];

			int count = 1;
			int uEmb = -1;
			for (int i = 0; i < cand.length; i++) {
				if (map[cand[i]] == 0) {
					if (rand.next() % count == 0) {
						uEmb = cand[i];
					}
					count++;
				}
			}
			if (uEmb < 0) {
				return false;
			}
			return addLeaf(map[vEmb], uEmb);
		}

		public boolean removeRandomLeaf() {
			if (leafs.size == 0) {
				return false;
			}
			int vEmb = leafs.array[rand.next() % leafs.size];
			return removeLeaf(vEmb);
		}

		public boolean addLeaf(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				return false;
			}

			lastOperation = OP_ADD;
			lastEmb = vEmb;
			lastV = v;

			map[vEmb] = v;
			leafs.add(vEmb);
			nodes.add(vEmb);
			embbedEdgeCount[vEmb]++;
			embbedEdgeCount[next]++;

			leafs.remove(next);
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}

		public boolean removeLeaf(int vEmb) {
			int v = map[vEmb];
			if (v == 0 || embbedEdgeCount[vEmb] != 1) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					System.out.println(u);
				}
				return false;
			}

			lastOperation = OP_REMOVE;
			lastEmb = vEmb;
			lastV = v;

			embbedEdgeCount[next]--;
			embbedEdgeCount[vEmb]--;

			leafs.remove(vEmb);
			nodes.remove(vEmb);
			map[vEmb] = 0;

			if (embbedEdgeCount[next] == 1) {
				leafs.add(next);
			} else if (embbedEdgeCount[next] == 0) {
				leafs.remove(next);
			}
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				dc[u][v]--;
				dc[v][u]--;
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount--;
				}
			}
			return true;
		}

		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
		
		public boolean isFull() {
			return edgeCount == E;
		}
	}
}


class Solver3 {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;

	private static final int ADD_RETRY = 10;

	private final Xor128 rand = new Xor128();
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;

	public Solver3(Data data) {

		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;

		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}

	public Result solve() {
		State state = new State();
		int len = Vemb / V;
		int L = (int)(Math.sqrt(Vemb) + 0.000001);
		Set<Integer> set = new HashSet<>();
		for (int i = 1; i <= V; i ++) {
			set.add(i);
		}

		for (int i = 0; i < V; i ++) {
			List<Integer> list = new ArrayList<>();
			for (int j = 0; j < len; j ++) {
				int idx = i * len + j;
				int h = (int)(Math.sqrt(idx) + 0.000001);
	
				int x = (idx) - h * h;
				int y = h;
				//i = 2, h=1, x=1, L = 5, y = 6
				if (x > h) {
					//i = 3, h = 1, L= 5, x = 2, y=
					y = h - (x - h);
					x = h;
				}
				if (h % 2 == 0) {
					int tmp = x;
					x = y;
					y = tmp;
				}
				int vEmb = y * L + x + 1;
				list.add(vEmb);
			}
			
			int max = -1;
			int maxV = 0;
			for (int v : set) {
				int now = 0;
				for (int vEmb : list) {
					for (int uEmb : h[vEmb]) {
						int u = state.map[uEmb];
						now += d[u][v] ? 1 : 0;
					}
				}
				if (now > max) {
					maxV = v;
					max = now;
				}
			}
			for (int vEmb : list) {
				state.addRoot(maxV, vEmb);
			}
			set.remove(maxV);
		}

		Result result = getResult(state, state.nodes.size, state.edgeCount);
		return result;
	}

	private Result getResult(State state, int... params) {
		Result result = new Result(state.map, state.score(), time, loop, params);
		return result;
	}

	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;

		// vEmb->v
		public final int[] map;

		private final IntSet leafs;
		private final IntSet nodes;
		private final int[] embbedEdgeCount;
		private final int[][] dc;

		private int edgeCount = 0;
		private int lastEmb = 0;
		private int lastV = 0;
		private int lastOperation = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			embbedEdgeCount = new int[Vemb + 1];
			leafs = new IntSet(Vemb + 1);
			nodes = new IntSet(Vemb + 1);
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			embbedEdgeCount = Arrays.copyOf(state.embbedEdgeCount, Vemb + 1);
			leafs = new IntSet(state.leafs);
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}

		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				addRoot(v, vEmb);
			}
		}

		public void addRoot(int v, int vEmb) {
			map[vEmb] = v;
			nodes.add(vEmb);

			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			
		}

		public boolean randomNext() {
			for (int k = 0; k < ADD_RETRY; k ++) {
				boolean ret = false;
				if (rand.next() % Xor128.RAND_INF < Xor128.RAND_INF) {
					ret = addRandomLeaf();
				} else {
					ret = removeRandomLeaf();
				}
				if (ret) {
					return true;
				}
			}
			return false;
		}

		public void undo() {
			if (lastOperation == OP_ADD) {
				removeLeaf(lastEmb);
			} else if (lastOperation == OP_REMOVE) {
				addLeaf(lastV, lastEmb);
			}
			lastOperation = 0;
		}

		public boolean addRandomLeaf() {
			if (nodes.size == 0) {
				return false;
			}
			int vEmb = nodes.array[rand.next() % nodes.size];
			int[] cand = h[vEmb];

			int count = 1;
			int uEmb = -1;
			for (int i = 0; i < cand.length; i++) {
				if (map[cand[i]] == 0) {
					if (rand.next() % count == 0) {
						uEmb = cand[i];
					}
					count++;
				}
			}
			if (uEmb < 0) {
				return false;
			}
			return addLeaf(map[vEmb], uEmb);
		}

		public boolean removeRandomLeaf() {
			if (leafs.size == 0) {
				return false;
			}
			int vEmb = leafs.array[rand.next() % leafs.size];
			return removeLeaf(vEmb);
		}

		public boolean addLeaf(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				return false;
			}

			lastOperation = OP_ADD;
			lastEmb = vEmb;
			lastV = v;

			map[vEmb] = v;
			leafs.add(vEmb);
			nodes.add(vEmb);
			embbedEdgeCount[vEmb]++;
			embbedEdgeCount[next]++;

			leafs.remove(next);
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}

		public boolean removeLeaf(int vEmb) {
			int v = map[vEmb];
			if (v == 0 || embbedEdgeCount[vEmb] != 1) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					System.out.println(u);
				}
				return false;
			}

			lastOperation = OP_REMOVE;
			lastEmb = vEmb;
			lastV = v;

			embbedEdgeCount[next]--;
			embbedEdgeCount[vEmb]--;

			leafs.remove(vEmb);
			nodes.remove(vEmb);
			map[vEmb] = 0;

			if (embbedEdgeCount[next] == 1) {
				leafs.add(next);
			} else if (embbedEdgeCount[next] == 0) {
				leafs.remove(next);
			}
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				dc[u][v]--;
				dc[v][u]--;
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount--;
				}
			}
			return true;
		}

		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
		
		public boolean isFull() {
			return edgeCount == E;
		}
	}
}

class Solver2 {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;

	private static final double START_TEMP = 100;
	private static final double END_TEMP = 10;
	private static final int ADD_RETRY = 10;

	private final Xor128 rand = new Xor128();
	private final long start = System.currentTimeMillis();
	private final long timeLimit;
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;

	public Solver2(Data data) {

		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;
		this.timeLimit = data.timeLimit;

		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}

	public Result solve() {
		State state = new State();
		State bestState = new State(state);
		int bestScore = bestState.score();
		int currentScore = state.score();

		int update = 0;
		int notUpdate = 0;
		int failed = 0;
		int epoch = 1;
		long start = 0;
		while (!isEnd()) {
			if (!state.randomNext()) {
				failed++;
				continue;
			}
			int nextScore = state.score();
			int scoreDiff = nextScore - currentScore;

			double temp = START_TEMP + (END_TEMP - START_TEMP) * (time - start) / (timeLimit - start);
			double probability = Math.exp((double) scoreDiff / temp);
			boolean forceNext = probability > (double) (rand.next() % Xor128.RAND_INF) / Xor128.RAND_INF;


			if (scoreDiff > 0 || forceNext) {
				currentScore = nextScore;
				update ++;
				if (bestScore < currentScore) {
					bestState = new State(state);
					bestScore = currentScore;
				}
			} else {
				notUpdate ++;
				if (notUpdate > 100000) {
					state = new State();
					notUpdate = 0;
					start = time - 1;
					epoch ++;
				}
				state.undo();
			}
		}
		Result result = getResult(bestState, state.nodes.size, bestState.edgeCount, epoch, update, failed);
		return result;
	}

	private Result getResult(State state, int... params) {
		Result result = new Result(state.map, state.score(), time, loop, params);
		return result;
	}

	private boolean isEnd() {
		loop++;
		if (loop % 100 != 0) {
			return false;
		}

		long end = System.currentTimeMillis();
		time = end - start;
		return time >= timeLimit;
	}

	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;

		// vEmb->v
		public final int[] map;

		private final IntSet leafs;
		private final IntSet nodes;
		private final int[] embbedEdgeCount;
		private final int[][] dc;

		private int edgeCount = 0;
		private int lastEmb = 0;
		private int lastV = 0;
		private int lastOperation = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			embbedEdgeCount = new int[Vemb + 1];
			leafs = new IntSet(Vemb + 1);
			nodes = new IntSet(Vemb + 1);
			init();
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			embbedEdgeCount = Arrays.copyOf(state.embbedEdgeCount, Vemb + 1);
			leafs = new IntSet(state.leafs);
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}

		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				map[vEmb] = v;
				nodes.add(vEmb);

				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					if (u == 0) {
						continue;
					}
					if (d[u][v] && dc[u][v] == 0) {
						edgeCount++;
					}
					dc[u][v]++;
					dc[v][u]++;
				}
			}
		}

		public boolean randomNext() {
			for (int k = 0; k < ADD_RETRY; k ++) {
				boolean ret = false;
				if (rand.next() % Xor128.RAND_INF < Xor128.RAND_INF) {
					ret = addRandomLeaf();
				} else {
					ret = removeRandomLeaf();
				}
				if (ret) {
					return true;
				}
			}
			return false;
		}

		public void undo() {
			if (lastOperation == OP_ADD) {
				removeLeaf(lastEmb);
			} else if (lastOperation == OP_REMOVE) {
				addLeaf(lastV, lastEmb);
			}
			lastOperation = 0;
		}

		public boolean addRandomLeaf() {
			if (nodes.size == 0) {
				return false;
			}
			int vEmb = nodes.array[rand.next() % nodes.size];
			int[] cand = h[vEmb];

			int count = 1;
			int uEmb = -1;
			for (int i = 0; i < cand.length; i++) {
				if (map[cand[i]] == 0) {
					if (rand.next() % count == 0) {
						uEmb = cand[i];
					}
					count++;
				}
			}
			if (uEmb < 0) {
				return false;
			}
			return addLeaf(map[vEmb], uEmb);
		}

		public boolean removeRandomLeaf() {
			if (leafs.size == 0) {
				return false;
			}
			int vEmb = leafs.array[rand.next() % leafs.size];
			return removeLeaf(vEmb);
		}

		public boolean addLeaf(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				return false;
			}

			lastOperation = OP_ADD;
			lastEmb = vEmb;
			lastV = v;

			map[vEmb] = v;
			leafs.add(vEmb);
			nodes.add(vEmb);
			embbedEdgeCount[vEmb]++;
			embbedEdgeCount[next]++;

			leafs.remove(next);
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}

		public boolean removeLeaf(int vEmb) {
			int v = map[vEmb];
			if (v == 0 || embbedEdgeCount[vEmb] != 1) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					System.out.println(u);
				}
				return false;
			}

			lastOperation = OP_REMOVE;
			lastEmb = vEmb;
			lastV = v;

			embbedEdgeCount[next]--;
			embbedEdgeCount[vEmb]--;

			leafs.remove(vEmb);
			nodes.remove(vEmb);
			map[vEmb] = 0;

			if (embbedEdgeCount[next] == 1) {
				leafs.add(next);
			} else if (embbedEdgeCount[next] == 0) {
				leafs.remove(next);
			}
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				dc[u][v]--;
				dc[v][u]--;
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount--;
				}
			}
			return true;
		}

		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
		
		public boolean isFull() {
			return edgeCount == E;
		}
	}
}

class Solver1 {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;

	private static final double START_TEMP = 1000;
	private static final double END_TEMP = 10;

	private final Xor128 rand = new Xor128();
	private final long start = System.currentTimeMillis();
	private final long timeLimit;
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;
	private double temp;

	public Solver1(Data data) {

		this.V = data.g.length - 1;
		this.Vemb = data.h.length - 1;
		this.g = data.g;
		this.h = data.h;
		this.timeLimit = data.timeLimit;

		d = new boolean[V + 1][V + 1];
		int e = 0;
		for (int v = 1; v <= V; v++) {
			e += g[v].length;
			for (int u : g[v]) {
				d[u][v] = d[v][u] = true;
			}
		}
		this.E = e;
	}

	public Result solve() {
		State state = new State();
		State bestState = new State(state);
		int bestScore = bestState.score();
		int currentScore = state.score();

		while (!isEnd()) {
			if (!state.randomNext()) {
				continue;
			}
			int nextScore = state.score();
			int scoreDiff = nextScore - currentScore;

			double probability = Math.exp(((double) scoreDiff) / temp);
			boolean forceNext = probability > (double) (rand.next() % Xor128.RAND_INF) / Xor128.RAND_INF;

			if (scoreDiff > 0 || forceNext) {
				currentScore = nextScore;
				if (bestScore < currentScore) {
					bestState = new State(state);
					bestScore = currentScore;
				}
			} else {
				state.undo();
			}
		}

		return getResult(bestState);
	}

	private Result getResult(State state) {
		Result result = new Result(state.map, state.score(), time, loop, new int[] {state.nodes.size});
		return result;
	}

	private boolean isEnd() {
		loop++;
		if (loop % 100 != 0) {
			return false;
		}

		long end = System.currentTimeMillis();
		time = end - start;
		temp = START_TEMP + (END_TEMP - START_TEMP) * time / timeLimit;
		return time >= timeLimit;
	}

	class State {
		public static final int OP_ADD = 1;
		public static final int OP_REMOVE = 2;

		// vEmb->v
		public final int[] map;

		private final IntSet leafs;
		private final IntSet nodes;
		private final int[] embbedEdgeCount;
		private final int[][] dc;

		private int edgeCount = 0;
		private int lastEmb = 0;
		private int lastV = 0;
		private int lastOperation = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			embbedEdgeCount = new int[Vemb + 1];
			leafs = new IntSet(Vemb + 1);
			nodes = new IntSet(Vemb + 1);
			init();
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			dc = new int[V + 1][V + 1];
			for (int i = 1; i <= V; i++) {
				for (int j = 1; j <= V; j++) {
					dc[i][j] = state.dc[i][j];
				}
			}
			embbedEdgeCount = Arrays.copyOf(state.embbedEdgeCount, Vemb + 1);
			leafs = new IntSet(state.leafs);
			nodes = new IntSet(state.nodes);
			edgeCount = state.edgeCount;
		}

		public void init() {
			for (int v = 1; v <= V; v++) {
				int vEmb;
				while (true) {
					vEmb = rand.next() % Vemb + 1;
					if (map[vEmb] == 0) {
						break;
					}
				}
				map[vEmb] = v;
				nodes.add(vEmb);

				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					if (u == 0) {
						continue;
					}
					if (d[u][v] && dc[u][v] == 0) {
						edgeCount++;
					}
					dc[u][v]++;
					dc[v][u]++;
				}
			}
		}

		public boolean randomNext() {
			boolean ret;
			if (rand.next() % Xor128.RAND_INF < Xor128.RAND_INF / 9 * 8) {
				ret = addRandomLeaf();
			} else {
				ret = removeRandomLeaf();
			}
			return ret;
		}

		public void undo() {
			if (lastOperation == OP_ADD) {
				removeLeaf(lastEmb);
			} else if (lastOperation == OP_REMOVE) {
				addLeaf(lastV, lastEmb);
			}
			lastOperation = 0;
		}

		public boolean addRandomLeaf() {
			if (nodes.size == 0) {
				return false;
			}
			int v = nodes.array[rand.next() % nodes.size];
			int[] cand = h[map[v]];
			int vEmb = cand[rand.next() % cand.length];
			return addLeaf(v, vEmb);
		}

		public boolean removeRandomLeaf() {
			if (leafs.size == 0) {
				return false;
			}
			int vEmb = leafs.array[rand.next() % leafs.size];
			return removeLeaf(vEmb);
		}

		public boolean addLeaf(int v, int vEmb) {
			if (map[vEmb] > 0) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				return false;
			}

			lastOperation = OP_ADD;
			lastEmb = vEmb;
			lastV = v;

			map[vEmb] = v;
			leafs.add(vEmb);
			nodes.add(vEmb);
			embbedEdgeCount[vEmb]++;
			embbedEdgeCount[next]++;

			leafs.remove(next);
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount++;
				}
				dc[u][v]++;
				dc[v][u]++;
			}
			return true;
		}

		public boolean removeLeaf(int vEmb) {
			int v = map[vEmb];
			if (v == 0 || embbedEdgeCount[vEmb] != 1) {
				return false;
			}
			int next = -1;
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == v) {
					if (next > 0) {
						return false;
					}
					next = uEmb;
				}
			}
			if (next < 0) {
				for (int uEmb : h[vEmb]) {
					int u = map[uEmb];
					System.out.println(u);
				}
				return false;
			}

			lastOperation = OP_REMOVE;
			lastEmb = vEmb;
			lastV = v;

			embbedEdgeCount[next]--;
			embbedEdgeCount[vEmb]--;

			leafs.remove(vEmb);
			nodes.remove(vEmb);
			map[vEmb] = 0;

			if (embbedEdgeCount[next] == 1) {
				leafs.add(next);
			} else if (embbedEdgeCount[next] == 0) {
				leafs.remove(next);
			}
			for (int uEmb : h[vEmb]) {
				int u = map[uEmb];
				if (u == 0) {
					continue;
				}
				dc[u][v]--;
				dc[v][u]--;
				if (d[u][v] && dc[u][v] == 0) {
					edgeCount--;
				}
			}
			return true;
		}

		public int score() {
			int score = POINT_BASE + edgeCount * POINT_EDGE - nodes.size * POINT_V + V;
			score += edgeCount == E ? POINT_FULL : 0;
			return score;
		}
	}
}
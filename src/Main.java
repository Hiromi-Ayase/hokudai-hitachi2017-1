import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

class Solver {

	private static final int POINT_BASE = 5000;
	private static final int POINT_EDGE = 100;
	private static final int POINT_FULL = 100000;
	private static final int POINT_V = 1;

	private static final double START_TEMP = 100;
	private static final double END_TEMP = 10;
	private static final int ADD_RETRY = 10;
	private static final double EPS = 0.000000001;
	private final Xor128 rand = new Xor128();
	private final long start = System.currentTimeMillis();
	private final long timeLimit;
	private final long timeLimit2 = 3000;
	private final int V, E, Vemb;
	private final int[][] g;
	private final int[][] h;
	private final boolean[][] d;
	private int loop = 0;
	private long time;
	private double temp;

	public Solver(Data data) {

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

	private static class FindPath {
		private final int[][] h;
		private final int n;

		public final int[] prev;
		public final int[] goal;
		public final int[] distance;

		public FindPath(int Vemb, int V, int[][] h) {
			this.h = h;
			this.n = Vemb + 1;
			this.prev = new int[n];
			this.distance = new int[n];
			this.goal = new int[V + 1];
		}

		public void dijk(int from, int[] map) {

			int I = Integer.MAX_VALUE / 2;
			Arrays.fill(distance, I);
			Arrays.fill(prev, -1);
			Arrays.fill(goal, -1);
			MinHeap q = new MinHeap(n);
			q.add(from, 0);
			distance[from] = 0;

			while (q.size() > 0) {
				int cur = q.argmin();
				q.remove(cur);

				for (int i = 0; i < h[cur].length; i++) {
					int next = h[cur][i];
					int nd = distance[cur] + 1;

					if (nd < distance[next]) {
						prev[next] = cur;
						distance[next] = nd;
						if (map[next] == 0) {
							q.update(next, nd);
						} else {
							goal[map[next]] = next;
						}
					}
				}
			}
		}
	}

	public Result solve() {

		FindPath fp = new FindPath(Vemb, V, h);
		int L = (int) (Math.sqrt(Vemb) + EPS);

		State state = solve2();
		IntSet now = new IntSet(Vemb + 1);

		State bestState = new State(state);
		int maxScore = bestState.score();

		while (!isEnd()) {
			loop++;
			int v = loop % V + 1;

			int min = Integer.MAX_VALUE;
			IntSet minSet = null;

			state.removeInv(v);

			for (int vEmb = 1; vEmb <= Vemb; vEmb++) {
				if (state.map[vEmb] > 0) {
					continue;
				}

				boolean flg = false;
				for (int uEmb : h[vEmb]) {
					if (state.map[uEmb] > 0) {
						flg = true;
						break;
					}
				}
				if (!flg)
					continue;

				now.clear();
				now.add(vEmb);

				int failed = 0;
				fp.dijk(vEmb, state.map);

				List<int[]> list = new ArrayList<>();
				for (int u : g[v]) {
					list.add(new int[] { u, fp.distance[u] });
				}
				Collections.sort(list, (o1, o2) -> o1[1] - o2[1]);

				for (int[] elem : list) {
					int u = elem[0];
					int p = fp.goal[u];
					if (p < 0) {
						failed++;
					} else {
						while (p > 0 && !now.contains(p)) {
							now.add(p);
							p = fp.prev[p];
						}
					}
				}

				int score = now.size + failed * L;
				if (score < min) {
					min = score;
					minSet = new IntSet(now);
				}
			}
			for (int i = 0; i < minSet.size; i++) {
				state.set(v, minSet.array[i]);
			}

			if (maxScore < state.score()) {
				maxScore = state.score();
				bestState = new State(state);
			}
		}
		// dump(state);
		Result result = getResult(bestState, bestState.nodes.size, bestState.edgeCount);
		return result;
	}

	public State solve2() {
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
		return bestState;
	}

	private State annealing(State state) {
		State bestState = new State(state);
		int bestScore = bestState.score();
		int currentScore = state.score();
		int L = (int) (Math.sqrt(Vemb) + 0.000001);
		while (!isEnd2()) {
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

	private boolean isEnd2() {
		loop++;
		if (loop % 100 != 0) {
			return false;
		}


		long end = System.currentTimeMillis();
		time = end - start;
		temp = START_TEMP + (END_TEMP - START_TEMP) * time / timeLimit2;
		return time >= timeLimit2;
	}


	private boolean isEnd() {
		loop++;
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
		private final int[] inv;
		private final int[][] dc;

		private int edgeCount = 0;

		public State() {
			map = new int[Vemb + 1];
			dc = new int[V + 1][V + 1];
			inv = new int[V + 1];
			nodes = new IntSet(Vemb + 1);
		}

		public State(State state) {
			map = Arrays.copyOf(state.map, Vemb + 1);
			inv = Arrays.copyOf(state.inv, V + 1);
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

		public boolean removeInv(int v) {
			int vEmb = inv[v];
			if (map[vEmb] != v) {
				return false;
			}

			ArrayDeque<Integer> queue = new ArrayDeque<>();
			queue.add(vEmb);
			remove(vEmb);
			while (queue.size() > 0) {
				vEmb = queue.poll();
				for (int uEmb : h[vEmb]) {
					if (map[uEmb] == v) {
						remove(uEmb);
						queue.add(uEmb);
					}
				}
			}
			return false;
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
			inv[v] = vEmb;
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

class MinHeap {
	public int[] a;
	public int[] map;
	public int[] imap;
	public int n;
	public int pos;
	public static int INF = Integer.MAX_VALUE;

	public MinHeap(int m) {
		n = m + 2;
		a = new int[n];
		map = new int[n];
		imap = new int[n];
		Arrays.fill(a, INF);
		Arrays.fill(map, -1);
		Arrays.fill(imap, -1);
		pos = 1;
	}

	public int add(int ind, int x) {
		int ret = imap[ind];
		if (imap[ind] < 0) {
			a[pos] = x;
			map[pos] = ind;
			imap[ind] = pos;
			pos++;
			up(pos - 1);
		}
		return ret != -1 ? a[ret] : x;
	}

	public int update(int ind, int x) {
		int ret = imap[ind];
		if (imap[ind] < 0) {
			a[pos] = x;
			map[pos] = ind;
			imap[ind] = pos;
			pos++;
			up(pos - 1);
		} else {
			a[ret] = x;
			up(ret);
			down(ret);
		}
		return x;
	}

	public int remove(int ind) {
		if (pos == 1)
			return INF;
		if (imap[ind] == -1)
			return INF;

		pos--;
		int rem = imap[ind];
		int ret = a[rem];
		map[rem] = map[pos];
		imap[map[pos]] = rem;
		imap[ind] = -1;
		a[rem] = a[pos];
		a[pos] = INF;
		map[pos] = -1;

		up(rem);
		down(rem);
		return ret;
	}

	public int min() {
		return a[1];
	}

	public int argmin() {
		return map[1];
	}

	public int size() {
		return pos - 1;
	}

	private void up(int cur) {
		for (int c = cur, p = c >>> 1; p >= 1 && a[p] > a[c]; c >>>= 1, p >>>= 1) {
			int d = a[p];
			a[p] = a[c];
			a[c] = d;
			int e = imap[map[p]];
			imap[map[p]] = imap[map[c]];
			imap[map[c]] = e;
			e = map[p];
			map[p] = map[c];
			map[c] = e;
		}
	}

	private void down(int cur) {
		for (int c = cur; 2 * c < pos;) {
			int b = a[2 * c] < a[2 * c + 1] ? 2 * c : 2 * c + 1;
			if (a[b] < a[c]) {
				int d = a[c];
				a[c] = a[b];
				a[b] = d;
				int e = imap[map[c]];
				imap[map[c]] = imap[map[b]];
				imap[map[b]] = e;
				e = map[c];
				map[c] = map[b];
				map[b] = e;
				c = b;
			} else {
				break;
			}
		}
	}
}

class Xor128 {
	public static final int RAND_INF = 10000000;
	private int x = 123456789, y = 362436069, z = 521288629, w = 88675123;

	public int next() {
		int t = (x ^ (x << 11));
		x = y;
		y = z;
		z = w;
		return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
	}

	public void shuffle(int[] array) {
		int n = array.length;
		for (int i = n - 1; i >= 1; i--) {
			int j = next() % i;
			int tmp = array[i];
			array[i] = array[j];
			array[j] = tmp;
		}
	}
}

class IntSet {
	public int[] array;
	private int[] inv;
	public int size = 0;

	public IntSet(IntSet obj) {
		int n = obj.array.length;
		array = Arrays.copyOf(obj.array, n);
		inv = Arrays.copyOf(obj.inv, n);
		size = obj.size;
	}

	public void clear() {
		Arrays.fill(inv, -1);
		size = 0;
	}

	public IntSet(int size) {
		array = new int[size];
		inv = new int[size];
		Arrays.fill(inv, -1);
	}

	public boolean contains(int value) {
		int p = inv[value];
		return (p >= 0 && p < size);
	}

	public void remove(int value) {
		int p = inv[value];
		if (p < 0 || p >= size)
			return;
		size--;
		int u = array[size];
		array[p] = u;
		inv[u] = p;
		inv[value] = -1;
	}

	public void add(int value) {
		if (inv[value] >= 0)
			return;
		array[size] = value;
		inv[value] = size;
		size++;
	}
}

class Result {
	public final int loop;
	public final long time;
	public final int score;
	public final int[] map;
	public final int[] params;

	public Result(int[] map, int score, long time, int loop, int[] params) {
		this.loop = loop;
		this.time = time;
		this.score = score;
		this.map = map;
		this.params = params;
	}
}

class Data {
	public final int[][] g;
	public final int[][] h;
	public final String type;
	public final long timeLimit;

	public Data(int[][] g, int[][] h, String type, long timeLimit) {
		this.type = type;
		this.g = g;
		this.h = h;
		this.timeLimit = timeLimit;
	}
}

public class Main {

	private static long TIME_LIMIT = 29000;

	public void run(long timeLimit) {
		int V = ni();
		int E = ni();
		int[] from = new int[E];
		int[] to = new int[E];

		for (int i = 0; i < E; i++) {
			from[i] = ni();
			to[i] = ni();
		}
		int Vemb = ni();
		int Eemb = ni();
		int[] fromEmb = new int[Eemb];
		int[] toEmb = new int[Eemb];
		for (int i = 0; i < Eemb; i++) {
			fromEmb[i] = ni();
			toEmb[i] = ni();
		}
		Data data = new Data(packU(V + 1, from, to), packU(Vemb + 1, fromEmb, toEmb), "Main", timeLimit);

		Result result = new Solver(data).solve();
		int[][] ret = new int[V + 1][Vemb + 1];
		int[] ptr = new int[V + 1];
		for (int i = 1; i <= Vemb; i++) {
			int v = result.map[i];
			ret[v][ptr[v]++] = i;
		}
		for (int i = 1; i <= V; i++) {
			out.print(ptr[i] + " ");
			for (int j = 0; j < ptr[i]; j++) {
				out.print(ret[i][j] + " ");
			}
			out.println();
		}
	}

	private static int[][] packU(int n, int[] from, int[] to) {
		int[][] g = new int[n][];
		int[] p = new int[n];
		for (int f : from)
			p[f]++;
		for (int t : to)
			p[t]++;
		for (int i = 0; i < n; i++)
			g[i] = new int[p[i]];
		for (int i = 0; i < from.length; i++) {
			g[from[i]][--p[from[i]]] = to[i];
			g[to[i]][--p[to[i]]] = from[i];
		}
		return g;
	}

	public static void main(String[] args) {
		System.err.println("Main mode.");
		reader = new java.io.BufferedReader(new java.io.InputStreamReader(System.in), 32768);
		new Main().run(TIME_LIMIT);
		out.flush();
	}

	private static java.io.PrintWriter out = new java.io.PrintWriter(System.out);
	private static java.util.StringTokenizer tokenizer = null;
	private static java.io.BufferedReader reader;

	public static String next() {
		while (tokenizer == null || !tokenizer.hasMoreTokens()) {
			try {
				tokenizer = new java.util.StringTokenizer(reader.readLine());
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		return tokenizer.nextToken();
	}

	private static int ni() {
		return Integer.parseInt(next());
	}
}
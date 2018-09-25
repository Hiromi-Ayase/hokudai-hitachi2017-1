
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Test {

	private final Xor128 rand = new Xor128();
	private static final Class<?> RECORD_CLASS = findSolver();

	private static final long INF = 10000;
	private static final int MAX_V = 500;
	private static final int MAX_E = 20000;
	private static final int MAX_L = 60;
	private static final int STEP = 30;
	private static final int EPOCH = 1;

	private static final int MULTI_THREAD = 4;
	private static final long TIME_LIMIT = 29000;
	private static final boolean RECORD = false;

	private String[][] bestData = new String[0][];
	private int bestDataPtr = 0;
	private final Class<?> solverClass;
	private final String recordFile;

	public static void main(String[] args) {
		if (RECORD) {
			new Test(findSolver()).record(EPOCH);
		} else {
			new Test(null).eval(EPOCH);
		}
	}

	private Test(Class<?> solverClass) {
		this.solverClass = solverClass;
		if (RECORD_CLASS == null) {
			this.recordFile = null;
		} else {
			this.recordFile = String.format("data/record_%s_%d.txt", RECORD_CLASS.getSimpleName(), TIME_LIMIT);
		}
	}

	private void eval(int n) {
		System.out.println("Run tests...");
		if (recordFile != null && Files.exists(Paths.get(recordFile))) {
			bestData = load();
		}
		long total = 0;
		long maxTotal = 0;
		long diffTotal = 0;
		for (int i = 1; i <= n; i++) {
			EpochResult ret = epoch(i);
			total += ret.total;
			maxTotal += ret.totalMax;
			diffTotal += ret.totalDiff;
		}
		String diffStr = diffTotal > Integer.MIN_VALUE ? String.format(", diff: %d", diffTotal / n) : "";
		System.out.printf("average: %d, max average: %d, ratio: %f%s%n", total / n, maxTotal / n,
				(double) total / maxTotal, diffStr);
	}

	private void record(int n) {
		System.out.println("Record tests...");
		try (BufferedWriter bw = Files.newBufferedWriter(Paths.get(recordFile))) {
			bw.write(TestResult.CSV_HEADER);
			bw.write("\n");
			for (int i = 0; i < n; i++) {
				EpochResult epochResult = epoch(i);
				for (TestResult result : epochResult.list) {
					bw.write(result.toCsv());
					bw.write("\n");
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private String[][] load() {
		List<String[]> list = new ArrayList<>();
		try (BufferedReader br = Files.newBufferedReader(Paths.get(recordFile))) {
			br.readLine();
			String line;
			while ((line = br.readLine()) != null) {
				list.add(line.split(","));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return list.stream().toArray(size -> new String[size][]);
	}

	private EpochResult epoch(int k) {
		List<List<Worker>> list = new ArrayList<>();
		for (int i = 1; i <= STEP; i += MULTI_THREAD) {
			List<Worker> currentWorkers = new ArrayList<>();
			for (int j = i; j < Math.min(i + MULTI_THREAD, STEP + 1); j++) {
				String[] best = bestDataPtr < bestData.length ? bestData[bestDataPtr++] : null;
				int[][][] g = buildRandomGraph();
				Data data = new Data(g[0], g[1], "Random", TIME_LIMIT);
				currentWorkers.add(new Worker(j, data, best, solverClass));
			}
			list.add(currentWorkers);
		}

		EpochResult epochResult = new EpochResult(k);
		try {
			ExecutorService service = Executors.newFixedThreadPool(MULTI_THREAD);
			TestResult lastResult = null;

			for (List<Worker> currentWorkers : list) {
				List<Future<TestResult>> futureResults = new ArrayList<>();
				int count = 0;
				for (Worker worker : currentWorkers) {
					count++;
					if (count == currentWorkers.size()) {
						lastResult = worker.call();
					} else {
						futureResults.add(service.submit(worker));
					}
				}
				for (Future<TestResult> future : futureResults) {
					epochResult.add(future.get());
					System.out.println(future.get());
				}
				epochResult.add(lastResult);
				System.out.println(lastResult);
			}
			service.shutdown();
		} catch (InterruptedException | ExecutionException e) {
			e.printStackTrace();
		}

		System.out.println(epochResult);
		;
		System.out.println();
		return epochResult;
	}

	private int[][] buildKingsGraph(int L) {
		List<int[]> list = new ArrayList<>();
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				int e = i * L + j;
				if (j < L - 1) {
					list.add(new int[] { e, e + 1 });
				}
				if (i < L - 1) {
					list.add(new int[] { e, e + L });
					if (j < L - 1) {
						list.add(new int[] { e, e + L + 1 });
					}
					if (j > 0) {
						list.add(new int[] { e, e + L - 1 });
					}
				}
			}
		}
		int n = list.size();
		int[] from = new int[n];
		int[] to = new int[n];
		int ptr = 0;
		for (int[] v : list) {
			from[ptr] = v[0] + 1;
			to[ptr] = v[1] + 1;
			ptr++;
		}
		return packU(L * L + 1, from, to);
	}

	private int[][][] buildRandomGraph() {
		int V = rand.next() % (MAX_V - 1) + 2;

		int E = rand.next() % (Math.min(V * (V - 1) / 2, MAX_E) - (V - 1)) + V - 1;
		assert (V - 1 <= E && E <= V * (V - 1) / 2);

		int[] func = new int[V];
		for (int i = 0; i < V; i++) {
			func[i] = i;
		}
		rand.shuffle(func);

		Set<Long> edges = new HashSet<>();

		for (int i = 1; i < V; i++) {
			int parent = rand.next() % i;
			int u = func[i];
			int v = func[parent];
			if (u > v) {
				int tmp = u;
				u = v;
				v = tmp;
			}
			edges.add(h(u, v));
		}

		for (int i = 0; i < E - (V - 1); i++) {
			while (true) {
				int u = rand.next() % (V - 1);
				int v = rand.next() % (V - u) + u;
				if (edges.contains(h(u, v)))
					continue;
				edges.add(h(u, v));

				break;
			}
		}

		int[] from = new int[E];
		int[] to = new int[E];
		int ptr = 0;
		for (long e : edges) {
			int[] edge = dh(e);
			from[ptr] = edge[0] + 1;
			to[ptr] = edge[1] + 1;
			ptr++;
		}

		int Vemb = 0;
		int L;
		while (true) {
			L = rand.next() % MAX_L + 1;
			Vemb = L * L;
			if (Vemb >= V) {
				break;
			}
		}
		int[][] g = packU(V + 1, from, to);
		int[][] kings = buildKingsGraph(L);
		return new int[][][] { g, kings };
	}

	private long h(int u, int v) {
		return u * INF + v;
	}

	private int[] dh(long h) {
		int v = (int) (h % INF);
		h /= INF;
		int u = (int) h;
		return new int[] { u, v };
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

	private static Class<?> findSolver() {
		Class<?> now = null;
		;
		for (int i = 1000; i >= 0; i--) {
			try {
				if ((now = Class.forName("Solver" + i)) != null) {
					break;
				}
			} catch (ClassNotFoundException e) {
				continue;
			}
		}
		return now;
	}
}

class Worker implements Callable<TestResult> {
	private final int step;
	private final Data data;
	private final String[] best;
	private final Class<?> solverClass;

	public Worker(int step, Data data, String[] best, Class<?> solverClass) {
		this.step = step;
		this.data = data;
		this.best = best;
		this.solverClass = solverClass;
	}

	@Override
	public TestResult call() {
		System.gc();
		Result result;
		if (solverClass == null) {
			result = new Solver(data).solve();
		} else {
			try {
				Method method = solverClass.getMethod("solve");
				Object instance = solverClass.getConstructor(Data.class).newInstance(data);
				result = (Result) method.invoke(instance);
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		TestResult testResult = new TestResult(step, result, data, best);
		return testResult;
	}
}

class TestResult {
	public static final String CSV_HEADER = "step,score,max,ratio,time,memory,V,E,Vemb,type,loop";
	public final int score;
	public final int resultScore;
	public final int max;
	public final long memory;
	public final long time;
	public final int V;
	public final int E;
	public final int Vemb;
	public final String type;
	public final int loop;
	public final int step;
	public final int diff;
	public final String params;

	public TestResult(int step, Result result, Data data, String[] bestData) {
		this.step = step;
		this.score = score(result, data);
		this.resultScore = result.score;
		this.memory = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / 1024 / 1024;
		this.time = result.time;
		this.V = data.g.length - 1;
		this.E = Arrays.stream(data.g).mapToInt(o -> o.length).sum() / 2;
		this.Vemb = data.h.length - 1;
		this.type = data.type;
		this.loop = result.loop;
		this.params = Arrays.toString(result.params);
		int max = 0;
		for (int i = 1; i <= V; i ++) {
			max += Math.min(8, data.g[i].length) * 100;
		}
		this.max = max;
		if (bestData != null) {
			this.diff = result.score - Integer.parseInt(bestData[1]);
		} else {
			this.diff = Integer.MIN_VALUE;
		}
	}

	@Override
	public String toString() {
		String diffStr = diff != Integer.MIN_VALUE ? String.format(", diff:%+6d", diff) : "";
		String diffScore = score != resultScore ? String.format("(!%d)", resultScore) : "";
		return String.format(
				"[%04d] score:%7d/%7d%s (%.3f)%s, time:%4dms, mem:%3dMB, V:%3d, E:%5d, Vemb:%4d, Type:%6s, loop:%8d, params:%s",
				step, score, max, diffScore, (double) score / max, diffStr, time, memory, V, E, Vemb, type, loop,
				params);
	}

	public String toCsv() {
		return String.format("%d,%d,%d,%f,%d,%d,%d,%d,%d,%s,%d", step, score, max, (double) score / max, time, memory,
				V, E, Vemb, type, loop);
	}

	private static int score(Result result, Data data) {
		int Vemb = data.h.length - 1;
		int V = data.g.length - 1;
		int[] map = result.map;

		int score = 5000;

		int[][] d = new int[V + 1][V + 1];
		for (int v = 1; v <= V; v++) {
			for (int u : data.g[v]) {
				d[u][v] = d[v][u] = 100;
			}
		}

		int count = 0;
		boolean[] used = new boolean[V + 1];
		for (int i = 1; i <= Vemb; i++) {
			int from = map[i];
			if (from == 0) {
				continue;
			}
			count++;
			boolean found = false;
			for (int next : data.h[i]) {
				int to = map[next];
				if (to == 0) {
					continue;
				}
				score += d[from][to];
				d[from][to] = d[to][from] = 0;
				if (from == to)
					found = true;
			}
			if (used[from] && !found) {
				return 0;
			}
			used[from] = true;
		}

		for (int i = 1; i <= V; i++) {
			if (!used[i]) {
				return 0;
			}
		}
		score += V - count;
		return score;
	}
}

class EpochResult {
	public final int epoch;
	public final List<TestResult> list;
	public long total = 0;;
	public long totalMax = 0;
	public long totalDiff = 0;

	public EpochResult(int epoch) {
		this.epoch = epoch;
		this.list = new ArrayList<>();
	}

	public void add(TestResult result) {
		this.list.add(result);
		this.total += result.score;
		this.totalMax += result.max;
		this.totalDiff += result.diff;
	}

	@Override
	public String toString() {
		String diffStr = totalDiff > Integer.MIN_VALUE ? String.format(", diff:%+5d", totalDiff) : "";
		return String.format("*** [Epoch %04d] score:%d/%d (%f)%s ***", epoch, total, totalMax,
				(double) total / totalMax, diffStr);

	}
}

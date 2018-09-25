import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Test2 {
	public static void main(String[] args) {
		List<int[]> list = new ArrayList<>();
		int r = 9;
		rec(0, r, -1, new int[r], list);
		for (int[] v : list) {
			System.out.println(Arrays.toString(v));
		}
//		for (int i = 0; i < 15; i ++) {
//			
//		}
//		
//		long[][] dp = new long[2][15 * 15 * 15 * 15];
//		
//		for (int i = 0; i < 4; i ++) {
//			int from = i % 2;
//			int to = (i + 1) % 2;
//			
//			
//		}
		
	}

	private static void rec(int depth, int maxDepth, int pre, int[] now, List<int[]> list) {
		
		if (depth == maxDepth) {
			list.add(Arrays.copyOf(now, maxDepth));
			return;
		}

		for (int i = 0; i < maxDepth; i ++) {
			if (i == pre) continue;
			now[depth] = i;
			rec(depth + 1, maxDepth, i, now, list);
		}
	}

	private static class Rand {
		private int x = 123456789, y = 362436069, z = 521288629, w = 88675123;

		public int next() {
			int t = (x ^ (x << 11));
			x = y;
			y = z;
			z = w;
			return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
		}

	}
}

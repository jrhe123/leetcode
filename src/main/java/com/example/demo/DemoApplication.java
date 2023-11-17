package com.example.demo;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class DemoApplication {
	
	
	public static void main(String[] args) {		
		Solution solution = new Solution();
	}
}

class Solution {
	
	
	
//	SESSION 2: find all duplicates in array, positive integer array
	public List<Integer> findDuplicates(int[] nums) {
		
	}
	
	
	
//	SESSION 3: number of distinct island
	public int numDistinctIsland(int[][] grid) {
		if (grid == null || grid.length == 0) {
			return 0;
		}
		
		Set<String> set = new HashSet<>();
		int m = grid.length;
		int n = grid[0].length;
		
		for(int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (grid[i][j] == 1) {
					String pathString = compuatePath(
							grid, i, j, m, n, "X"
							);
					set.add(pathString);
				}
			}
		}
		
		
		return set.size();
	}
	
	
	// X = start
	// O = out of bound OR water
	// U = up
	// D = down
	// L = left
	// R = right
	private String compuatePath(
			int[][] grid, int i, int j, int m, int n, String sign
			) {
		// exit condition
		if (i < 0 || j < 0 || i >= m || j >= n || grid[i][j] == 0)
			return "O";
		
		grid[i][j] = 0;
		// recursive
		String leftString = compuatePath(
				grid, i, j - 1, m, n, "L"
				);
		String upString = compuatePath(
				grid, i - 1, j, m, n, "U"
				);
		String rightString = compuatePath(
				grid, i, j + 1, m, n, "R"
				);
		String downString = compuatePath(
				grid, i + 1, j, m, n, "D"
				);
		
		return leftString + upString + rightString + downString;
	}
	
	
	
//	SESSION 4: a/e/i/o/u
	public String reverseVowelString(String s) {
		char[] arr = s.toCharArray();
		
		int left = 0;
		int right = arr.length - 1;
		
		while (left < right) {
			boolean leftIsVowel = isVowel(arr[left]);
			boolean rightIsVowel = isVowel(arr[right]);
			
			if (leftIsVowel && rightIsVowel) {
				swapChar(
						arr, left, right
						);
				++left;
				--right;
			}
			
			if (!leftIsVowel) {
				++left;
			}
			if (!rightIsVowel) {
				--right;
			}
		}
		
		return new String(arr);
	}
	
	private boolean isVowel(char c) {
		char lowerCase = Character.toLowerCase(c);
		return lowerCase == 'a' || lowerCase == 'e' || lowerCase == 'i' ||
				lowerCase == 'o' || lowerCase == 'u';
	}
	
	private void swapChar(char[] arr, int i, int j) {
		char temp = arr[i];
		arr[i] = arr[j];
		arr[j] = temp;
	}
	
	
	
//	SESSION 5: reverse integer
	public int reserverInteger(int x) {
		long result = 0;
		
		while(x != 0) {
			
			int remainder = x % 10;
			result = result * 10 + remainder;
			
			// 2^31
			if (result > Integer.MAX_VALUE ||
					result < Integer.MIN_VALUE) {
				return 0;
			}
			
			x /= 10;
		}
		
		return (int) result;
	}
	
	
	
//	SESSION 6: maximum area of island
	public int maxAreaOfIsland(int[][] grid) {
		int max = 0;
		int n = grid.length;
		int m = grid[0].length;
		
		for(int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (grid[i][j] == 1) {
					int area = getArea(
							grid, i, j, n, m
							);
					max = Math.max(max, area);
				}
			}
		}
		
		return max;
	}
	
	private int getArea(int[][] grid, int i, int j, int n, int m) {
		// exit condition
		if (i < 0 || j < 0 || i >= n || j >= m || grid[i][j] == 0) {
			return 0;
		}
		// set it to zero
		grid[i][j] = 0;
		// recursive
		int left = getArea(
				grid, i, j - 1, n, m
				);
		int up = getArea(
				grid, i - 1, j, n, m
				);
		int right = getArea(
				grid, i, j + 1, n, m
				);
		int down = getArea(
				grid, i + 1, j, n, m
				);
		return left + up + right + down + 1;
	}
	
	
//	SESSION 7: first not repeat char
	public int firstUniqueCharacter(String s) {
		// 1. hashmap store key-value pair
		Map<Character, Integer> countMap = new HashMap<>();
		for(int i = 0; i < s.length(); i++) {
			
			char c = s.charAt(i);
			if (countMap.containsKey(c)) {
				countMap.put(c, countMap.get(c) + 1);
			} else {
				countMap.put(c, 1);
			}
		}
		// 2. iterator the hashmap
		for(int i = 0; i < s.length(); i++) {
			char c = s.charAt(i);
			
			if (countMap.get(c) == 1) {
				return i;
			}
		}
		return -1;
	}
}
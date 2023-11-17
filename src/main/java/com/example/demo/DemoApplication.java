package com.example.demo;

import java.util.ArrayList;
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
	
//	SESSION 1: dynamic programming, only move right / down in a 2-D array
	public int uniquePaths(int m, int n) {
		int[][] dp = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (i == 0 || j == 0) {
					dp[i][j] = 1;
				} else {
					dp[i][j] = dp[i-1][j] + dp[i][j-1];
				}
			}
		}
		
		return dp[m - 1][n-1];
	}
	
	
	
	
//	SESSION 2: find all duplicates in array, positive integer array
	public List<Integer> findDuplicates(int[] nums) {
		List<Integer> resultIntegers = new ArrayList<>();
		
		for (int i = 0; i < nums.length; i++) {
			int index = Math.abs(nums[i]) - 1;
			
			if (nums[index] < 0) {
				resultIntegers.add(index + 1);
			}
			
			nums[index] = -nums[index];
		}
		
		return resultIntegers;
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
	
	
	
//	SESSION 8: valid anagram - same characters in different order
//	ASCII table: a -> 97, b -> 98
	public boolean isValidAnagram(String s, String t) {
		int[] letters = new int[26];
		char[] arr1 = s.toCharArray();
		char[] arr2 = t.toCharArray();
		
		for (int i = 0; i < arr1.length; i++) {
			letters[arr1[i] - 'a']++;
		}
		for (int i = 0; i < arr2.length; i++) {
			letters[arr2[i] - 'a']--;
		}
		
		for (int i = 0; i < letters.length; i++) {
			if (letters[i] != 0) return false;
		}
		
		return true;
	}
	
	
//	SESSION 9: number of islands
	public int numberOfIsland(int[][] grid) {
		
		if (grid == null || grid.length == 0) return 0;
		
		int islandCount = 0;
		int m = grid.length;
		int n = grid[0].length;
		
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (grid[i][j] == 1) {
					traverseIsland(
							grid, i, j, m, n
							);
					++islandCount;
				}
			}
		}		
		
		return islandCount;
	}
	
	private void traverseIsland(
			int[][] grid, int i, int j, int m, int n
			) {
		// exit condition
		if (i < 0 || j < 0 || i >= m || j >= n || grid[i][j] == 0) return;
		
		grid[i][j] = 0;
		// recursive
		// left
		traverseIsland(
				grid, i, j-1, m, n
				);
		// up
		traverseIsland(
				grid, i-1, j, m, n
				);
		// right
		traverseIsland(
				grid, i, j+1, m, n
				);
		// down
		traverseIsland(
				grid, i + 1, j, m, n
				);
	}
	
	
}
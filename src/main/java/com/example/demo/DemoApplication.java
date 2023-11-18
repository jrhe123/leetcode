package com.example.demo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

import org.apache.el.util.ConcurrentCache;

import com.fasterxml.jackson.annotation.JacksonInject.Value;

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
	
	
//	SESSION 10: sum of two integers to target, array is sorted
	public int[] twoSum(int[] numbers, int target) {
		int left = 0;
		int right = numbers.length - 1;
		
		while(left < right) {
			int sum = numbers[left] + numbers[right];
			
			if (sum == target) return new int[] {
					left, right
			};
			else if (sum > target) --right;
			else ++left;
		}
		
		return null;
	}
	
	
//	SESSION 11: clone graph
	public Node cloneGraph(Node node) {
		if (node == null) return null;
		
		Map<Integer, Node> map = new HashMap<>();
		
		return cloneGraph(node, map);
	}
	
	//	DFS
	private Node cloneGraph(Node node, Map<Integer, Node> map) {
		
		if (map.containsKey(node.val)) return map.get(node.val);
		
		Node copyNode = new Node(node.val);
		map.put(node.val, copyNode);
		
		for (Node neighborNode : node.neighbors) {
			copyNode.neighbors.add(
					cloneGraph(neighborNode, map)
					);
		}
		
		return copyNode;
	}
	

//	SESSION 12: Sum of left leaves: all leaves in a given binary tree
	public int sumOfLeftLeaves(TreeNode root) {
		if(root == null) return 0;
		
		int sum = 0;
		if (root.leftTreeNode != null) {
			if (isLeaf(root.leftTreeNode)) sum += root.leftTreeNode.val;
			else {
				sum += sumOfLeftLeaves(root.leftTreeNode);
			}
		}
		
		sum += sumOfLeftLeaves(root.rightTreeNode);
		
		return sum;
	}
	
	private boolean isLeaf (TreeNode node) {
		return node.leftTreeNode == null && node.rightTreeNode == null;
	}
	
	
//	QUEUE
	public int sumOfLeftLeavesV2(TreeNode root) {
		if(root == null) return 0;
		
		int sum = 0;
		
		Queue<TreeNode> queue = new LinkedList<>();
		queue.add(root);
		
		while (!queue.isEmpty()) {
			TreeNode node = queue.poll();
			
			if(node.leftTreeNode != null) {
				if (isLeaf(node)) {
					sum += node.leftTreeNode.val;
				} else {
					queue.add(node.leftTreeNode);
				}
			}
			
			if (node.rightTreeNode != null)
				queue.add(node.rightTreeNode);
		}
		
		return sum;
	}
	
	
	
//	SESSION 13: Min stack, find minimum value in the stack
//	MinStack
	
	
	
//	SESSION 14: Find first missing positive integer, required O(1)
//	int test[] = {7,8,2000};
//	int findFirstMissingPositiveIntegerV1 = solution.findFirstMissingPositiveIntegerV1(test);
//	System.out.println(findFirstMissingPositiveIntegerV1);
	
//	O(N)
	public int findFirstMissingPositiveIntegerV1(int[] arr) {
		Set<Integer> set = new HashSet<>();
		for (int i = 0; i < arr.length; i++) {
			set.add(arr[i]);
		}
		for (int i = 0; i < arr.length + 1; i++) {
			if (!set.contains(i + 1)) {
				return i + 1;
			}
		}
		return arr.length;
	}
	
//	O(1)
	public int findFirstMissingPositiveIntegerV2(int[] nums) {
		if (nums == null || nums.length == 0) return 1;
		
		int n = nums.length, containsOne = 0;
		
		// step 1
		for (int i = 0; i < n; i++) {
			if (nums[i] == 1) {
				containsOne = 1;
			} else if (nums[i] <= 0 || nums[i] > n) {
				nums[i] = 1;
			}
		}
		if (containsOne == 0) return 1;
		
		// step 2
		for (int i = 0; i < n; i++) {
			int index = Math.abs(nums[i]) - 1;
			if (nums[index] > 0) nums[index] *= -1;
		}
		
		// step 3
		for (int i = 0; i < n; i++) {
			if (nums[i] > 0) {
				return i + 1;
			}
		}
		return n + 1;
	}
	
	
//	SESSION 15: Rotate image: 90 degree
	public void rotate(int[][] matrix) {
		if (matrix == null || matrix.length == 0 ) return;
		
		int n = matrix.length;
		transposeMatrix(matrix, n);
		reverseMatrix(matrix, n);
	}
	
	
	private void transposeMatrix(int[][] matrix, int n) {
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				swap(matrix, i, j, j, i);
			}
		}
	}
	
	private void reverseMatrix(int[][] matrix, int n) {
		for (int i = 0; i < n; i++) {
			int left = 0;
			int right = n - 1;
			while (left < right) {
				swap(matrix, i, left, i, right);
				
				++left;
				--right;
			}
		}
	}
	
	private void swap(int[][] matrix, int i, int j, int k, int l) {
		int temp = matrix[i][j];
		matrix[i][j] = matrix[k][l];
		matrix[k][l] = temp;
	}

	
	
//	SESSION 16: Lowest Common Ancestor of a binary tree
//	Post-order: from bottom to up
//	L -> R -> Visit Node
	private TreeNode resulTreeNode;
	
	public TreeNode lowestCommonAncestor(TreeNode root, TreeNode p, TreeNode q) {
		findLCA(root, p, q);
		return resulTreeNode;
	}
	
	private boolean findLCA(TreeNode root, TreeNode p, TreeNode q) {
		// exit condition
		if (root == null) return false;
		
		// recursive
		boolean left = findLCA(root.leftTreeNode, p, q);
		boolean right = findLCA(root.rightTreeNode, p, q);	
		boolean curr = root == p || root == q;
		
		if((left && right) || (left && curr) || (right && curr)) {
			resulTreeNode = root;
		}
		return left || right || curr;
	}
	
	
	
//	SESSION 17: Add Strings: "2859" + "293"
//	ASCII: subtract char '0'
	public String addStrins(String num1, String num2) {
		int i = num1.length() - 1;
		int j = num2.length() - 1;
		
		int carry = 0;
		StringBuilder resultBuilder = new StringBuilder();
		
		while (i > -1 || j > -1) {
			 int d1 = i > -1 ? num1.charAt(i) - '0' : 0;
			 int d2 = j > -1 ? num1.charAt(j) - '0' : 0;
			 
			 int sum = d1 + d2 + carry;
			 resultBuilder.append(sum % 10);
			 carry = sum / 10;
			 
			 --i;
			 --j;
		}
		
		if (carry == 1) resultBuilder.append(1);
		
		return resultBuilder.reverse().toString();
	}
	
	
//	SESSION 18: Min time difference
//	["03:00", "01:00", "23:30"]
//	180, 60, 1410
//	1440 sec a day
//	bucket sort
	public int findMinDifference(List<String> timePoints) {
		
		boolean[] bucket = new boolean[1440];
		
		for (String timePoint: timePoints) {
			String[] tStrings = timePoint.split(":");
			int hours = Integer.parseInt(tStrings[0]);
			int minutes = Integer.parseInt(tStrings[1]);
			
			int total = hours * 60 + minutes;
			
			if (bucket[total]) return 0;
			
			bucket[total] = true;
		}
		
		int min = Integer.MAX_VALUE;
		int first = -1;
		int prev = -1;
		int curr = -1;
		
		for (int i = 0; i < bucket.length; i++) {
			if (bucket[i]) {
				if (prev == -1) {
					prev = i;
					first = i;
				} else {
					curr = i;
					min = Math.min(min, curr - prev);
					prev = curr;
				}
			}
		}
		
		// check back circle to the second day
		int last = 1440 - curr + first;
		
		return Math.min(min, last);
	}
	
	
//	SESSION 19: Top K Frequent words
//	["i", "love", "leetcode", "i", "love", "coding"]
//	time: nlogn
//	space: n
//	Priority Queue
	public List<String> topKFrequent(String[] words, int k) {
		Map<String, Integer> map = new HashMap<>();
		
		// record down the frequency of each word
		for (String word: words) {
			map.put(word, map.getOrDefault(word, 0) + 1);
		}
		
		// prepare priority queue
		PriorityQueue<String> pQueue = new PriorityQueue<>(new Comparator<String>() {
			@Override
			public int compare(String s1, String s2) {
				// sort by frequency
				int frequence1 = map.get(s1);
				int frequence2 = map.get(s2);
				
				// sort by alphabetic
				if (frequence1 == frequence2) {
					return s2.compareTo(s1);
				}
				return frequence1 - frequence2;
			}
		});
		
		for (Map.Entry<String, Integer> entyEntry: map.entrySet()) {
			pQueue.add(entyEntry.getKey());
			
			// remove the less frequencies ones based on k
			if (pQueue.size() > k) pQueue.poll();
		}
		
		List<String> resultList = new ArrayList<>();
		while (!pQueue.isEmpty()) {
			resultList.add(pQueue.poll());
		}
		
		// reverse it back
		Collections.reverse(resultList);
		
		return resultList;
	}
	
	
//	SESSION 20: Merge two sorted lists
//  L1: 1 -> 2 -> 4 -> 10
//  L2: 1 -> 3 -> 4
	public ListNode mergeTwoLists(ListNode l1, ListNode l2) {
		
		// dummy head
		ListNode headListNode = new ListNode(-1);
		ListNode currListNode = headListNode;
		
		while (l1 != null || l2 != null) {
			
			if (l1 == null) {
				
				currListNode.nextListNode = l2;
				l2 = l2.nextListNode;
				
			} else if (l2 == null) {
				
				currListNode.nextListNode = l1;
				l1 = l1.nextListNode;
				
			} else if (l1.val < l2.val) {
				
				currListNode.nextListNode = l1;
				l1 = l1.nextListNode;
				
			} else { // l2.val < l1.val
				
				currListNode.nextListNode = l2;
				l2 = l2.nextListNode;
				
			}
			
			currListNode = currListNode.nextListNode;
		}
		
		// remove the dummy one
		return headListNode.nextListNode;
	}
	
	
//	SESSION 21: Making a large island
	private int[][] directions = {
			{1,0}, {-1, 0},
			{0,1}, {0, -1},
			};
	
	public int largestIsland(int[][] grid) {
		if (grid == null || grid.length == 0) return 0;
		
		int max = 0;
		int isIslandId = 2;
		int m = grid.length;
		int n = grid[0].length;
		
		Map<Integer, Integer> map = new HashMap<>();
		// label the existing island with different ID, and store it into map
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if(grid[i][j] == 1) {
					int size = getIslandSize(grid, i, j, isIslandId);
					
					max = Math.max(max, size);
					map.put(isIslandId, size);
					++isIslandId;
				}
			}
		}
		// test the new island with maximum
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (grid[i][j] == 0) {
					Set<Integer> set = new HashSet<>();
					for (int[] direction: directions) {
						// get the new position
						int x = direction[0] + i ;
						int y = direction[1] + j ;
						
						if (x > -1 && y > -1 &&
								x < m && y < n &&
								grid[x][y] != 0) {
							set.add(grid[x][y]);
						}
					}
					
					int sum = 1;
					for (int num : set) {
						int value = map.get(num);
						sum += value;
					}
					
					max = Math.max(sum, max);
				}
			}
		}
		
		return max;
	}
	
	private int getIslandSize(int[][] grid, int i, int j, int isIslandId) {
		// exit condition
		if (i < 0 || j < 0 ||
				i >= grid.length || j >= grid[0].length ||
				grid[i][j] != 1) return 0;
		
		grid[i][j] = isIslandId;
		
		// recursive
		int left = getIslandSize(grid, i, j - 1, isIslandId);
		int up = getIslandSize(grid, i - 1, j, isIslandId);
		int right = getIslandSize(grid, i, j + 1, isIslandId);
		int down = getIslandSize(grid, i + 1, j, isIslandId);
		
		return left + up + right + down + 1;
	}
	
	
//	SESSION 22: Text Justification
//	Greedy
//	
//	maxWidth = 15
//	words = ["What", "must", "be", "acknowledgment", "shall", "be"];
//	
//	output:
//	What___must__be (middle justify, left spaces > right spaces)
//	acknowledgment_ (left justify - single word)
//	shall_be_______ (left justify - last line)
	public List<String> fullJustify(String[] words, int maxWidth) {
		List<String> resultList = new ArrayList<>();
		
		int i = 0;
		int n = words.length;
		while (i < n) {
			
			int j = i + 1;
			int lineLength = words[i].length();
			
			while (j < n && 
					(lineLength + words[j].length() + (j - i - 1)) < maxWidth
					) {
				
				lineLength +=  words[j].length();
				++j;
			}
			
			int diff = maxWidth - lineLength;
			int numberOfWords = j - i;
			
			// check it's left / middle justify
			if (numberOfWords == 1 || j >= n) {
				resultList.add(
						leftJustify(words, diff, i, j)
						);
			} else {
				resultList.add(
						middleJustify(words, diff, i, j)
						);
			}
			
			// next line
			i = j;			
		}
		return resultList;
	}
	
	private String leftJustify(String[] words, int diff, int i, int j) {
		int spacesOnRight = diff - (j - i - 1);
		
		StringBuilder resuBuilder = new StringBuilder(words[i]);
		for (int k = i + 1; k < j; ++k) {
			resuBuilder.append(" " + words[k]);
		}
		
		resuBuilder.append(" ".repeat(spacesOnRight));
		
		return resuBuilder.toString();
	}
	
	private String middleJustify(String[] words, int diff, int i, int j) {
		int spacesNeeded = j - i - 1;
		
		int spaces = diff / spacesNeeded;
		int extraSpaces = diff % spacesNeeded;
		
		StringBuilder resuBuilder = new StringBuilder(words[i]);
		for (int k = i + 1; k < j; ++k) {
			extraSpaces--;
			int spacesToApply = spaces + (extraSpaces > 0 ? 1 : 0);
			
			resuBuilder.append(" ".repeat(spacesToApply) + words[k]);
		}
		return resuBuilder.toString();
	}
	
	
//	SESSION 23: Minimum Path Sum | Dynamic Programming
//  only right / down
//	2D grid
//	1,3,1
//	1,5,1
//	4,2,1
	public int minPathSum(int[][] grid) {
		int n = grid.length;
		int m = grid[0].length;
		
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < m; j++) {
				int top = i - 1 < 0 ? Integer.MAX_VALUE : grid[i-1][j];
				int left = j - 1 < 0 ? Integer.MAX_VALUE : grid[i][j-1];
				
				int min = top == Integer.MAX_VALUE && left == Integer.MAX_VALUE ?
						0 : Math.min(top, left);
				grid[i][j] += min;
			}
		}
		return grid[n - 1][m - 1];
	}
	
	
//	SESSION 24: Word Ladder | Breadth First Search
//	BFS
//	start: be
//	end:   ko
//	words: ["ce", "mo", "ko", "me", "co"]
	public int ladderLength(String beginWord, String endWord, List<String> wordList) {
		
		Set<String> set = new HashSet<>(wordList);
		if (!set.contains(endWord)) return 0;
		
		Queue<String> queue = new LinkedList<>();
		queue.add(beginWord);
		
		Set<String> visitedSet = new HashSet<>();
		visitedSet.add(beginWord);
		
		int changes = 1;
		
		while (!queue.isEmpty()) {
			
			int size = queue.size();
			for(int i = 0; i < size; i++) {
				String wordString = queue.poll();
				
				// exit condition
				if (wordString.equals(endWord)) return changes;
								
				for (int j = 0; j < wordString.length(); j++) {
					for (int k = 'a'; k < 'z'; k++) {
						char[] arr = wordString.toCharArray();
						arr[j] = (char) k;
						
						String string = new String(arr);
						if (set.contains(string) && 
								!visitedSet.contains(string)) {
							
							queue.add(string);
							visitedSet.add(string);
						}
					}
				}
			}
			++changes;
		}
		
		return 0;
	}
	
	
	
//	SESSION 25: Container With Most Water
//	heights = [1,8,6,2,5,4,8,3,7]
//	area: (R-L) * min(height[L], height[R])
//	move the shorter index of R or L
//	if same hegit, move both R & L
//	until L == R
	public int maxArea(int[] heights) {
		int n = heights.length;
		int left = 0;
		int right = n - 1;
		int max = 0;
		
		while (left < right) {
			int lh = heights[left];
			int rh = heights[right];
			
			int width = right - left;
			int height = Math.min(lh, rh);
			
			int area = width * height;
			max = Math.max(max, area);
			
			if (lh < rh) {
				left++;
			} else if (lh > rh) {
				right--;
			} else {
				left++;
				right--;
			}
		}
		
		return max;
	}
	
//	SESSION 26: Union Find Data Structure 
//	Number of Connected Components in an Undirected Graph
//	
//	Connect the possible disjoint subset, and then
//	find how many components in total
//	n = 5
//	[0,1,2,3,4]
//	edges = [[0,1], [1,2]. [3,4]]
//	
//	index: 		[0,1,2,3,4]
//	=> parents: [1,2,2,4,4]
//	components: when index == parent
//	=> parents: [2,2,2,4,4]
	public int countComponents(int n, int[][] edges) {
		int[] ids = new int[n];
		
		// initial: pointing to themselves as parent
		for (int i = 0; i < ids.length; i++) {
			ids[i] = i;
		}
		
		// union the edges
		for (int[] edge: edges) {
			union(edge[0], edge[1], ids);
		}
		
		// find the parent of edge
		Set<Integer> set = new HashSet<>();
		for (int i = 0; i < ids.length; i++) {
			// we use set for duplicate elimination
			set.add(
					find(i, ids)
					);
		}
		
		// number of components
		return set.size();
	}
	
	// find the parent of edge1, edge2, and point one of parent to other one
	private void union(int edge1, int edge2, int[] ids) {
		int parent1 = find(edge1, ids);
		int parent2 = find(edge2, ids);

		ids[parent1] = parent2;
	}
	
	// return parent of the edge
	private int find(int edge, int[] ids) {
		if (ids[edge] != edge) {
			// not parent
			// recursive: find the parent
			ids[edge] = find(ids[edge], ids);
		}
		return ids[edge];
	}
	
	
//	SESSION 27: Sliding Window Algorithm
//	Longest Substring Without Repeating Characters
//	
//	PWWKEW
//	=> WKE, KEW => 3
//	index i & j, starting at left
	public int lengthOfLongestSubstring(String s) {
		if (s == null || s.length() == 0) return 0;
		
		int i = 0;
		int j = 0;
		int max = 0;
		
		Set<Character> set = new HashSet<>();
		while (i < s.length()) {
			
			char c = s.charAt(i);
			
			while (set.contains(c)) {
				set.remove(s.charAt(j));
				++j;
			}
			
			set.add(c);
			max = Math.max(max, i - j + 1);
			++i;
		}
		
		return max;
	}
	

//	SESSION 28: Trie Data Structure
//	Tree
//	insert / search / startsWith
//	Trie

	

//	SESSION 29: Reorder Data in log files
//	letter log / digital log
//	letter log: let1 art can
//	digit log: dig1 8 1 5 1
	
//	split the first space: 
//	
//	letter logs come before digit log
//	letter logs sorted lexicographically, not including the id
//	letter logs are equal, sort by id lexicographically
//	digit log maintain their order
	public String[] recordLogFiles(String[] logs) {
		Arrays.sort(logs, (log1, log2) -> {
			// log1 < log2 => -1
			// log1 > log2 => 1
			// log1 == log2 => 0
			
			int index1 = log1.indexOf(" ");
			String id1 = log1.substring(0, index1);
			String main1 = log1.substring(index1 + 1);
			
			int index2 = log2.indexOf(" ");
			String id2 = log2.substring(0, index2);
			String main2 = log2.substring(index2 + 1);
			
			boolean isDigit1 = Character.isDigit(main1.charAt(0));
			boolean isDigit2 = Character.isDigit(main2.charAt(0));
			
			if (!isDigit1 && !isDigit2) {
				int value = main1.compareTo(main2);
				if (value == 0) {
					return id1.compareTo(id2);
				}
				return value;
			}
			
			return isDigit1 ? (isDigit2 ? 0 : 1) : -1;		
		});
		
		return logs;
	}
	
	
//	SESSION 30: Integer to Roman Numerals
//	I   V   X   L   C   D   M
//	1	5	10	50	100	500	1000
//	
//	edge case: IV: 4
//	IV	IX	XL	XC	CD	 CM
//	4	9	40	90	400	 900
	
//	range: 1 ~ 3999
	
	
	
	
	
	
	
	
	
}


class Trie {
	
	private Node root;
	
	public Trie () {
		root = new Node('\0');
	}
	
	public void insert (String word) {
		Node currNode = root;
		for (int i = 0; i < word.length(); i++) {
			char c = word.charAt(i);
			if(currNode.children[c - 'a'] == null) {
				currNode.children[c - 'a'] = new Node(c);
			}
			// move down
			currNode = currNode.children[c - 'a'];
		}
		// mark the last c - isWord = true
		currNode.isWord = true;
	}
	
	public boolean search(String word) {
		Node node = getNode(word);
		return node != null && node.isWord;
	}
	
	public boolean startsWith(String prefix) {
		Node node = getNode(prefix);
		return node != null;
	}
	
	private Node getNode(String word) {
		Node curr = root;
		for (int i = 0; i < word.length(); i++) {
			char c = word.charAt(i);
			if (curr.children[c - 'a'] == null) return null;
			// move down
			curr = curr.children[c - 'a'];
		}
		// found it
		return curr;
	}
	
	class Node {
		public char c;
		public boolean isWord;
		public Node[] children;
		
		public Node(char _c) {
			this.c = _c;
			isWord = false;
			children = new Node[26];
		}
		
	}
}



class ListNode {
	int val;
	ListNode nextListNode;
	
	ListNode() {}
	
	ListNode(int _val) {
		this.val  = _val;
	}
	
	ListNode(int _val, ListNode _nextListNode) {
		this.nextListNode  = _nextListNode;
	}
	
}


class MinStack {
	
	private Stack<Integer> s;
	private Stack<Integer> t;
	
	public MinStack() {
		s = new Stack<>();
		t = new Stack<>();
	}
	
	public void push(int x) {
		s.push(x);
		int min  = t.isEmpty() || x < t.peek() ? x : t.peek();
		t.push(min);
	}
	
	public void pop() {
		s.pop();
		t.pop();
	}
	
	public int top() {
		return s.peek();
	}
	
	public int getMin() {
		return t.peek();
	}
}


// binary tree node
class TreeNode {
	int val;
	TreeNode leftTreeNode;
	TreeNode rightTreeNode;
	
	TreeNode(int x) {
		val = x;
	}
}


class Node {
	public Node() {
		val = 0;
		neighbors = new ArrayList<>();
	}
	
	public Node(int _val) {
		val = _val;
		neighbors = new ArrayList<>();
	}
	
	public Node(int _val, ArrayList<Node> _neighbors) {
		val = _val;
		neighbors = _neighbors;
	}
	
	public int val;
	public List<Node> neighbors;
}

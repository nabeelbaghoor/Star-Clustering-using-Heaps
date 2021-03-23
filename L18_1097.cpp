#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include <math.h>
#include <list>
using namespace std;
struct Point {
	double x;
	double y;
	Point(double a, double b)
		:x(a), y(b)
	{}
	bool operator<(Point &obj)
	{
		return (x<obj.x);
	}
	bool operator>(Point& obj)
	{
		return (x > obj.x);
	}
};
struct Pair {
	int a;
	int b;
	Pair(int _a, int _b)
		:a(_a), b(_b)
	{}
	void swap(int& a, int& b)
	{
		swap(a, b);
	}
	Pair()
	{
	}
};
double distance(Point A, Point B) {
	return sqrt(pow((B.x - A.x), 2) + pow((B.y - A.y), 2));
}

template <typename k, typename v>
struct HeapItem
{
	k key;
	v value;
	HeapItem(k _key, v _value)
		:key(_key), value(_value)
	{
	}
	void swap(HeapItem<k, v>& a, HeapItem<k, v>& b)
	{
		swap(a.key, b.key);
		swap(a.value, b.value);
	}
	HeapItem()
	{
	}
};
template <typename k, typename v>
class MinHeap
{
	HeapItem<k, v>* arr;
	int capacity;
	int totalItems;
public:
	MinHeap()
	{
		capacity = 0;
		arr = nullptr;
		totalItems = 0;
	}
	MinHeap(int _capacity)
	{
		capacity = _capacity;
		arr = new HeapItem <k, v>[capacity];
		totalItems = 0;
	}
	void Allocate(int _capacity)
	{
		capacity = _capacity;
		arr = new HeapItem <k, v>[capacity];
		totalItems = 0;
	}
	void insert(k _key, v _value)
	{
		if (totalItems == capacity)
		{
			HeapItem<k, v>* temp = new HeapItem<k, v>[capacity * 2];
			for (int i = 0; i < capacity; i++)
				temp[i] = arr[i];
			delete arr;
			arr = temp;
		}
		arr[++totalItems].key = _key;
		arr[totalItems].value = _value;
		int i = totalItems;
		while ((i > 1) && (arr[i].key < arr[i / 2].key))
		{
			swap(arr[i], arr[i / 2]);
			i = i / 2;
		}
	}
	void getMin(v& _value)//gives value
	{
		if (!isEmpty())
			_value = arr[1].value;
	}
	bool isEmpty() const
	{
		return totalItems == 0;
	}
	bool isFull() const
	{
		return totalItems ==capacity;
	}
	void Heapify(int i)
	{
		int max = i;
		if (2 * i <= totalItems && arr[2 * i].key < arr[max].key)
			max = 2 * i;
		if (2 * i + 1 <= totalItems && arr[2 * i + 1].key < arr[max].key)
			max = 2 * i + 1;
		if (i != max)
		{
			swap(arr[i], arr[max]);
			Heapify(max);

		}
	}
	void deleteMin()
	{
		if (!isEmpty())
		{
			arr[1] = arr[totalItems];
			totalItems--;
			Heapify(1);
		}
	}
	void deleteMin(k& _key,v& _value)
	{
		if (!isEmpty())
		{
			_key = arr[1].key;
			_value = arr[1].value;
			arr[1] = arr[totalItems];
			totalItems--;
			Heapify(1);
		}
	}
	void decreaseKey(int index, k _key, v _value)//decrease key with value
	{
		if (index > 1 && index <= totalItems && _key < arr[index].key)
		{
			arr[index].key = _key;
			arr[index].value = _value;
			while (index > 1 && arr[index].key < arr[index / 2].key)
			{
				swap(arr[index], arr[index / 2]);
				index = index / 2;
			}
		}
	}
	void decreaseKey(int index,k _key)//only key
	{
		if (index > 1 && index <= totalItems && _key < arr[index].key)
		{
			arr[index].key = _key;
			while (index>1 && arr[index].key<arr[index/2].key)
			{
				swap(arr[index], arr[index / 2]);
				index = index / 2;
			}
		}
	}
	void Remove(int index)//key or value
	{
		if (index > 1 && index <= totalItems)
		{
			decreaseKey(index, arr[1].key-1, arr[1].value);//LESS KEY ONLY
			deleteMin();
		}
	}
	HeapItem<k,v>&operator[](int i)
	{
		return arr[i];
	}
	int Size()
	{
		return totalItems;
	}
	void BuildHeap(int _totalItems)
	{
		totalItems = _totalItems-1;
		for (int i = totalItems/2; i>=1; i--)
			Heapify(i);
	}
	void BuildHeap()
	{
		for (int i = totalItems / 2; i >= 1; i--)
			Heapify(i);
	}	
};
void PrintHeap(MinHeap<double, Pair>& heap)
{
	cout<<"Heap:\n";
	for (int i = 1; i <= heap.Size(); i++)
		cout << i << "\t" << heap[i].key << "\t(" << heap[i].value.a << "," << heap[i].value.b << ")\n";
}
void BuildMatrix(MinHeap<double, Pair>& heap, int** matrix)
{
	for (int m = 1; m <=heap.Size(); m++)
	{
		matrix[heap[m].value.a][heap[m].value.b] = m;
		matrix[heap[m].value.b][heap[m].value.a] = m;
	}
}
void PrintMatrixIndex(MinHeap<double, Pair>& heap, int** matrix,int size)
{
	cout << "Matrix with Index:\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << matrix[i][j]<<"\t";
		cout << endl;
	}
}
void PrintMatrixDistance(MinHeap<double, Pair>& heap, int** matrix, int size)
{
	cout << "Matrix with Distance:\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == 0 || j == 0)
				cout << matrix[i][j] << "\t";
			else if (matrix[i][j] != -1)
				cout << heap[matrix[i][j]].key << "\t";
			else
				cout << -1 << "\t";
		}
		cout << endl;
	}
}
void read(string filename)
{
	ifstream fin;
	fin.open(filename);
	if (fin.is_open())
	{
		double x, y;
		int N, M;
		char junk;
		fin >> N;
		fin >> junk;
		fin >> M;
		vector<list<Point>> groups(N + 1);

		//create 2D matrix
		int** matrix;
		matrix = new int* [N + 1];
		for (int i = 0; i < N + 1; i++)
			matrix[i] = new int[N + 1]();
		//diagonal Entries already set to zero,as all entries zero

		//to make matrix readable
		for (int i = 0; i < N+1; i++)
		{
			matrix[0][i] = i;
			matrix[i][0] = i;
		}
		int k = 1;
		while (!fin.eof()) {
			fin >> x;
			fin >> y;
			groups[k++].push_back(Point(x, y));
		}
		int heapSize = ((N * N) - N) / 2;
		MinHeap<double, Pair>heap(heapSize+1);//
		heap[0].key = -2.0;
		//FillHeap
		k = 1;
		for (int m = 1; m < N + 1; m++)
			for (int n = m + 1; n < N + 1; n++)
					heap[k++] = HeapItem<double, Pair>(distance(groups[m].front(), groups[n].front()), Pair(m, n));
		
		//buildHeap
		cout << "heapSize: " << heapSize << endl;
		heap.BuildHeap(k);
		
		//buildMatrix
		BuildMatrix(heap,matrix);

		double dist=0.0;
		double MinDist;
		Pair p;
		int noOfGroups = N;
		int w = 1;
		PrintMatrixIndex(heap, matrix, N + 1);
		//PrintMatrixDistance(heap, matrix, N + 1);
		vector<bool>flags(N + 1);
		for (int i = 0; i < N + 1; i++)
			flags[i] = 1;
		cout << "Heap Before Extracting Min:\n";
		PrintHeap(heap);
		while (!heap.isEmpty() && noOfGroups>M)
		{
			heap.deleteMin(dist,p);
			if(flags[p.a] && flags[p.b])
			{ 
				flags[p.b] = 0;
				BuildMatrix(heap,matrix);//must build it
				cout << "Heap Extracted: " << p.a << "," << p.b << "\t" << "dist: " << dist << endl;
				groups[p.a].merge(groups[p.b]);
				//MergeColumns
				for (int i = 1; i <= p.a; i++)
				{
					if (matrix[i][p.a] != 0 && matrix[i][p.b] != 0 && matrix[i][p.a] != -1 && matrix[i][p.b]!=-1)
					{
						MinDist = ((heap[matrix[i][p.a]].key < heap[matrix[i][p.b]].key) ? heap[matrix[i][p.a]].key : heap[matrix[i][p.b]].key);
						heap[matrix[i][p.a]] = HeapItem<double, Pair>(MinDist, Pair(i, p.a));
						heap.BuildHeap();
						heap.Remove(matrix[i][p.b]);

						if (matrix[i][p.b])//if Not Zero
							matrix[i][p.b] = -1;
						BuildMatrix(heap,matrix);
					}
				}
				//MergeRows
				for (int i = p.a + 1; i < N + 1; i++)
				{
					if (matrix[p.a][i] != 0 && matrix[p.b][i] != 0 && matrix[p.a][i] != -1 && matrix[p.b][i]!=-1)
					{
						MinDist = ((heap[matrix[p.a][i]].key < heap[matrix[p.b][i]].key) ? heap[matrix[p.a][i]].key : heap[matrix[p.b][i]].key);
						heap[matrix[p.a][i]] = HeapItem<double, Pair>(MinDist, Pair(p.a,i));
						heap.BuildHeap();
						heap.Remove(matrix[p.b][i]);

						if(matrix[p.b][i])//if Not Zero
							matrix[p.b][i] = -1;
						BuildMatrix(heap,matrix);
					}
				}
				for (int i = 1; i < N + 1; i++)//for columns
					if(matrix[i][p.b])
						matrix[i][p.b] = -1;
				for (int i = 1; i < N + 1; i++)//for rows
					if (matrix[p.b][i])
						matrix[p.b][i] = -1;

				noOfGroups--;
				cout << "ITTERATION#" << w++ << endl;
				PrintMatrixIndex(heap, matrix, N + 1);
				//PrintMatrixDistance(heap, matrix, N + 1);
				
				if (!heap.isEmpty() && noOfGroups > M)
				{
					cout << "Heap Before Extracting Min:\n";
					PrintHeap(heap);
				}

			}
		}
		cout << "\nOUTPUT:\n";
		int  num = 1;
		for (int i = 1; i < groups.size(); i++)
		{
			if (!groups[i].empty())
			{
				cout << "Group" << num++<<": ";
				for (list<Point>::iterator itr = groups[i].begin(); itr != groups[i].end(); itr++)
					cout << "(" << (*itr).x << "," << (*itr).y << ")"<<" ";
				cout << endl;
			}
		}
	}
}
int main()
{
	read("a.txt");
	system("pause");
}
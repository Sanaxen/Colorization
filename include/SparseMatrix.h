#ifndef _SparseMatrix_H
#define _SparseMatrix_H

#include <vector>
#include <utility>
#include <map>
#include <unordered_map>

#if 10


class SparseMatrix
{
	int N;
	std::map<std::pair<int, int>, int> sparseMatrix;

public:
	SparseMatrix(const int n)
	{
		N = n;
	}

	inline void set( int val,  int row,  int col)
	{
		sparseMatrix[std::pair<int, int>(row, col)] = val;
	}
	inline int get(int row, int col) const
	{
		const std::pair<int, int> key = std::pair<int, int>(row, col);
		
		std::map<std::pair<int, int>, int>::const_iterator ret = sparseMatrix.find(key);

		if (ret != sparseMatrix.end())
		{
			return ret->second;
		}
		return -1;
	}

};

#else
class SparseMatrix
{
	int m, n;

	std::vector<int> * vals, *cols, *rows;

public:
	SparseMatrix(const int n)
	{
		this->construct(n, n);
	}
	~SparseMatrix(void)
	{
		if (this->vals != NULL) {
			delete this->vals;
			delete this->cols;
		}

		delete this->rows;
	}
	void construct(const int rows, const int columns)
	{
		if (rows < 1 || columns < 1) {
			throw "Matrix dimensions cannot be zero or negative.";
		}

		this->m = rows;
		this->n = columns;

		this->vals = this->cols = NULL;
		this->rows = new std::vector<int>(rows + 1, 1);
	}
	int get(const int row, const int col) const
	{

		const int x = this->rows->at(row - 1) - 1;
		const int n = this->rows->at(row) - 1;

		for (int i = x; i < n; i++) {

			int actual = this->cols->at(i);

			if (actual == col) {
				return this->vals->at(i);

			}
			else if (actual > col) {
				break;
			}
		}

		return 0;
	}
	SparseMatrix & set(const int val, const int row, const int col)
	{
		int pos = this->rows->at(row - 1) - 1;
		int actual = -1;

		const int n = this->rows->at(row) - 1;
		for (; pos < n; pos++) {
			actual = this->cols->at(pos);

			if (actual == col) {
				break;

			}
			else if (actual > col) {
				break;
			}
		}
		if (actual != col) {
			if (val != 0) {
				this->insert(pos, row, col, val);
			}

		}
		else if (val == 0) {
			this->remove(pos, row);
		}
		else {
			this->vals->at(pos) = val;
		}
		return *this;
	}

	void insert(const int index, const int row, const int col, int val)
	{
		if (this->vals == NULL) {
			this->vals = new std::vector<int>(1, val);
			this->cols = new std::vector<int>(1, col);

		}
		else {
			this->vals->insert(this->vals->begin() + index, val);
			this->cols->insert(this->cols->begin() + index, col);
		}

		const int mm = m;
		for (int i = row; i <= m; i++) {
			this->rows->at(i) = this->rows->at(i) + 1;
		}
	}


	void remove(const int index, const int row)
	{
		this->vals->erase(this->vals->begin() + index);
		this->cols->erase(this->cols->begin() + index);

		const int mm = m;
		for (int i = row; i <= mm; i++) {
			this->rows->at(i) = this->rows->at(i) - 1;
		}
	}
};
#endif

#endif
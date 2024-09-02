/* Copyright Andrei-Marian Rusanescu 311CAb 2023-2024 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

/* constant for % big results on various mathematical operations */
#define kt 10007

/* Allocates memory in the initial stage */
void allocate_memory(int ****container, int **rows, int **columns, int *dim)
{
	*dim = 1;
	*container = (int ***)malloc(*dim * sizeof(int **));
	if (!(*container)) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	*rows = (int *)malloc(*dim * sizeof(int));
	if (!(*rows)) {
		free(*container);
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	*columns = (int *)malloc(*dim * sizeof(int));
	if (!(*columns)) {
		free(*container);
		free(*rows);
		fprintf(stderr, "malloc() failed\n");
		return;
	}
}

/* Reallocates if the matrix array is fully loaded or half empty */
void reallocate_memory(int ****container, int **rows, int **columns, int dim)
{
	*container = (int ***)realloc(*container, dim * sizeof(int **));
	*rows = (int *)realloc(*rows, dim * sizeof(int));
	*columns = (int *)realloc(*columns, dim * sizeof(int));

	if (!(*container) || !(*rows) || !(*columns)) {
		fprintf(stderr, "realloc() failed\n");
		return;
	}
}

/* Allocates memory for square matrices */
void all_mat(int ***mat, int n)
{
	int **a = (int **)malloc(n * sizeof(int *));
	if (!a) {
		fprintf(stderr, "malloc() failed\n");
		*mat = NULL;
		return;
	}

	*mat = a;

	for (int i = 0; i < n; i++) {
		a[i] = (int *)malloc(n * sizeof(int));
		if (!a[i]) {
			fprintf(stderr, "malloc() for line %d failed\n", i);
			for (int j = i - 1; j >= 0; j++)
				free(a[j]);
			free(a);
			*mat = NULL;
			return;
		}
	}
}

/* Frees the allocated memory for square matrices */
void free_mat(int **mat, int n)
{
	for (int i = 0; i < n; i++)
		free(mat[i]);
	free(mat);
}

/* takes advantage of pointers' abilities, swaps two matrices */
void swap_matrices(int ***a, int ***b)
{
	int **tmp = *a;
	*a = *b;
	*b = tmp;
}

/* Swaps the matrices' indexes if swap_matrices is used */
void swap_pos(int *ln1, int *ln2)
{
	int tmp = *ln1;
	*ln1 = *ln2;
	*ln2 = tmp;
}

void check_positive(int *a)
{
	if (*a < 0)
		*a += kt;
}

/* Logarithmic complexity exponentiation algorithm */
void power_matrix(int ***mat, int n, int k)
{
	/* returns the identity matrix if power is 0 */
	if (k == 0) {
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				if (i == j)
					(*mat)[i][j] = 1;
				else
					(*mat)[i][j] = 0;
			}
		return;
	}
	int **a, **b, **c;
	all_mat(&a, n), all_mat(&b, n), all_mat(&c, n);
	if (!a || !b || !c) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j)
				b[i][j] = 1;
			else
				b[i][j] = 0;
			a[i][j] = 0, c[i][j] = 0;
		}
	}
	while (k > 0) {// power is an odd number
		if (k % 2 == 1) {
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					for (int l = 0; l < n; ++l) {
						a[i][j] += (b[i][l] * (*mat)[l][j]) % kt;
						check_positive(&a[i][j]);
					} //every time power is odd, it multiplies by mat: 1 * mat
					a[i][j] %= kt;
					check_positive(&a[i][j]);
				}
			// a is used to temporarily store the values, then emptied
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j) {
					b[i][j] = a[i][j];
					a[i][j] = 0;
				}
		} // if power is an even number
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				for (int l = 0; l < n; ++l) {// mat * mat
					c[i][j] += ((*mat)[i][l] * (*mat)[l][j]) % kt;
					check_positive(&c[i][j]);
				}
				c[i][j] %= kt;
				check_positive(&c[i][j]);
			} // mat becomes mat ^ 2, then mat ^ 4 and so on
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				(*mat)[i][j] = c[i][j];
				c[i][j] = 0;
			} // x^n = x^(n/2) * x^(n/2);
		k /= 2;
	} // the logarithmic exponentiation of matrices is done in place
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j)
			(*mat)[i][j] = b[i][j];
	} // frees all resources allocated;
	free_mat(a, n), free_mat(b, n), free_mat(c, n);
}

/* Function to calculate sum of the matrix' elements */
void find_sum(int ***container, int **s, int cnt, int *ln, int *cn)
{
	for (int k = 0; k < cnt; ++k) {
		for (int i = 0; i < ln[k]; ++i) {
			for (int j = 0; j < cn[k]; ++j) {
				(*s)[k] += (container)[k][i][j] % kt;
				check_positive(&(*s)[k]);
			}
		}
		(*s)[k] %= kt;
		check_positive(&(*s)[k]);
	}
}

/* Loads matrix in the matrix array */
void read_matrix(int ****container, int *rows, int *columns, int *cnt)
{
	int m, n;
	scanf("%d %d", &m, &n);
	rows[*cnt] = m;
	columns[*cnt] = n;

	(*container)[*cnt] = (int **)malloc(m * sizeof(int *));
	if (!(*container)[*cnt]) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}

	for (int i = 0; i < m; ++i) {
		(*container)[*cnt][i] = (int *)malloc(n * sizeof(int));
		if (!(*container)[*cnt][i]) {
			fprintf(stderr, "malloc() failed\n");
			return;
		}
		for (int j = 0; j < n; ++j)
			scanf("%d", &(*container)[*cnt][i][j]);
	}
	(*cnt)++;
}

/* Prints matrix sizes */
void command_D(int cnt, int *ln, int *cn)
{
	int index;
	scanf("%d", &index);
	if (cnt < index + 1 || index < 0)
		printf("No matrix with the given index\n");
	else
		printf("%d %d\n", ln[index], cn[index]);
}

/* Prints a specific matrix */
void command_P(int ***container, int *ln, int *cn, int cnt)
{
	int index;
	scanf("%d", &index);
	if (cnt < index + 1 || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		for (int i = 0; i < ln[index]; ++i) {
			for (int j = 0; j < cn[index]; ++j)
				printf("%d ", (container)[index][i][j]);
			printf("\n");
		}
	}
}

/* Concatenates matrix in place, using specific indications */
void command_C(int ****container, int *ln, int *cn, int cnt)
{
	int index;
	scanf("%d", &index);
	if (cnt < index + 1 || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		int rows_C, columns_C;// rows and columns for the new matrix
		scanf("%d", &rows_C);
		int *li = (int *)malloc(rows_C * sizeof(int));
		if (!li) {
			fprintf(stderr, "malloc() failed\n");
			return;
		} // which rows new matrix gets from old matrix;
		for (int i = 0; i < rows_C; ++i)
			scanf("%d", &li[i]);

		scanf("%d", &columns_C);
		int *co = (int *)malloc(columns_C * sizeof(int));
		if (!co) {
			fprintf(stderr, "malloc() failed\n");
			return;
		} // which columns new matrix gets from old matrix;
		for (int i = 0; i < columns_C; ++i)
			scanf("%d", &co[i]);

		// auxiliary matrix wherein the changes will first be made
		int **aux = (int **)malloc(rows_C * sizeof(int *));
		if (!aux) {
			fprintf(stderr, "malloc() failed\n");
			return;
		}

		for (int i = 0; i < rows_C; ++i) {
			aux[i] = (int *)malloc(columns_C * sizeof(int));
			if (!aux[i]) {
				fprintf(stderr, "malloc() failed\n");
				for (int j = i - 1; j >= 0; --j)
					free(aux[j]);
				free(aux);
				return;
			}
		}
		// aux has the given sizes;

		for (int i = 0; i < rows_C; ++i)
			for (int j = 0; j < columns_C; ++j)
				aux[i][j] = (*container)[index][li[i]][co[j]];
		// swapping the old matrix and sizes with the new ones;
		swap_matrices(&(*container)[index], &aux);
		swap_pos(&ln[index], &rows_C);
		swap_pos(&cn[index], &columns_C);

		for (int i = 0; i < rows_C; ++i)
			free(aux[i]);
		free(aux), free(li), free(co);
	}
}

/* Multiplies two matrices */
void command_M(int ****cont, int **rows, int **columns, int *cnt, int *nr)
{
	int i1, i2;
	scanf("%d %d", &i1, &i2);
	if (*cnt < i1 + 1 || *cnt < i2 + 1 || i1 < 0 || i2 < 0) {
		printf("No matrix with the given index\n");
	} else {
		if ((*columns)[i1] != (*rows)[i2]) {
			printf("Cannot perform matrix multiplication\n");
		} else {
			int **a = (int **)malloc((*rows)[i1] * sizeof(int *));
			if (!a) {
				fprintf(stderr, "malloc() failed\n");
				return;
			}
			for (int i = 0; i < (*rows)[i1]; ++i) {
				a[i] = (int *)calloc((*columns)[i2], sizeof(int));
				if (!a[i]) {
					fprintf(stderr, "calloc() failed\n");
					for (int j = i - 1; j >= 0; --j)
						free(a[j]);
					free(a);
					return;
				}
			}
			// matrix multiplication:
			for (int i = 0; i < (*rows)[i1]; ++i)
				for (int k = 0; k < (*columns)[i2]; ++k) {
					for (int j = 0; j < (*rows)[i2]; ++j) {
						a[i][k] += ((*cont)[i1][i][j] * (*cont)[i2][j][k]) % kt;
						check_positive(&a[i][k]);
					}
					a[i][k] %= kt;
					check_positive(&a[i][k]);
				}
			// the matrix is now made, needs to be loaded up in the array;
			// checks if needed to double the size of the array;
			if (*cnt == *nr) {
				(*nr) *= 2;
				reallocate_memory(cont, rows, columns, *nr);
			}
			// updates dimensions of the newly created matrix;
			(*rows)[*cnt] = (*rows)[i1];
			(*columns)[*cnt] = (*columns)[i2];

			// loading it up
			(*cont)[*cnt] = (int **)malloc((*rows)[i1] * sizeof(int *));
			if (!(*cont)[*cnt]) {
				fprintf(stderr, "malloc() failed\n");
				return;
			}

			for (int i = 0; i < (*rows)[i1]; ++i) {
				(*cont)[*cnt][i] = (int *)malloc((*columns)[i2] * sizeof(int));
				if (!(*cont)[*cnt][i]) {
					fprintf(stderr, "malloc() failed\n");
					return;
				}
				for (int j = 0; j < (*columns)[i2]; ++j)
					(*cont)[*cnt][i][j] = a[i][j];
			}
			for (int i = 0; i < (*rows)[i1]; ++i)
				free(a[i]);
			free(a);
			(*cnt)++;
		}
	}
}

/* Orders matrices in an ascending order,
 * by the sum of the matrix elements */
void command_O(int ****container, int *rows, int *columns, int cnt)
{
	int *s = (int *)calloc((cnt), sizeof(int));
	if (!s) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	find_sum(*container, &s, cnt, rows, columns);
	// bubble sort;
	for (int i = 0; i < cnt - 1; ++i) {
		for (int j = 0; j < cnt - i - 1; ++j) {
			if (s[j] > s[j + 1]) {
				swap_pos(&s[j], &s[j + 1]);
				swap_matrices(&(*container)[j], &(*container)[j + 1]);
				swap_pos(&rows[j], &rows[j + 1]);
				swap_pos(&columns[j], &columns[j + 1]);
			}
		}
	}
	free(s);
}

/* Transposes matrix in place */
void command_T(int ****container, int **rows, int **columns, int cnt)
{
	int index;
	scanf("%d", &index);
	if (cnt < index + 1 || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		int rowsT = (*columns)[index], columnsT = (*rows)[index];
		int **t = (int **)malloc(rowsT * sizeof(int *));
		if (!t) {
			fprintf(stderr, "malloc() failed\n");
			return;
		}
		for (int i = 0; i < rowsT; ++i) {
			t[i] = (int *)malloc(columnsT * sizeof(int));
			if (!t[i]) {
				fprintf(stderr, "malloc() failed\n");
				for (int j = i - 1; j >= 0; --j)
					free(t[j]);
				free(t);
			} //f orms the transpose
			for (int j = 0; j < columnsT; ++j)
				t[i][j] = (*container)[index][j][i];
		}
		swap_matrices(&(*container)[index], &t);
		swap_pos(&(*rows)[index], &rowsT);
		swap_pos(&(*columns)[index], &columnsT);
		for (int i = 0; i < rowsT; ++i)
			free(t[i]);
		free(t);
	}
}

/* Raises a matrix to a given power */
void command_R(int ****container, int *rows, int *columns, int cnt)
{
	int index;
	scanf("%d", &index);
	if (cnt < index + 1 || index < 0) {
		printf("No matrix with the given index\n");
	} else {
		int power;
		scanf("%d", &power);
		if (power < 0) {
			printf("Power should be positive\n");
		} else {
			if (rows[index] != columns[index]) {
				printf("Cannot perform matrix multiplication\n");
			} else {
				// logarithmic exponentiation;
				power_matrix(&(*container)[index], rows[index], power);
			}
		}
	}
}

/* Sums two square matrices of the same size */
void add(int **a, int **b, int size, int ***s)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j) {
			(*s)[i][j] = (a[i][j] + b[i][j]) % kt;
			check_positive(&(*s)[i][j]);
		}
}

/* Subtracts two square matrices of the same size */
void subtract(int **a, int **b, int size, int ***s)
{
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j) {
			(*s)[i][j] = (a[i][j] - b[i][j]) % kt;
			check_positive(&(*s)[i][j]);
		}
}

/* Strassen algorithm for matrix multiplication; - recursive */
void strs(int **mat1, int **mat2, int sz, int ***prod)
{
	if (sz == 1) {
		(*prod)[0][0] = (mat1[0][0] * mat2[0][0]) % kt;
		check_positive(&(*prod)[0][0]);
		return; // default condition for recursion;
	} // the matrix a(i) and b(i) are used in the multiplication process;
	int sz2 = sz / 2, **a1, **a2, **a3, **a4, **b1, **b2, **b3, **b4;
	int **m1, **m2, **m3, **m4, **m5, **m6, **m7;
	all_mat(&a1, sz2); all_mat(&a2, sz2); all_mat(&a3, sz2); all_mat(&a4, sz2);
	all_mat(&b1, sz2); all_mat(&b2, sz2); all_mat(&b3, sz2); all_mat(&b4, sz2);
	all_mat(&m1, sz2); all_mat(&m2, sz2); all_mat(&m3, sz2); all_mat(&m4, sz2);
	all_mat(&m5, sz2); all_mat(&m6, sz2); all_mat(&m7, sz2);
	if (!a1 || !a2 || !a3 || !a4 || !b1 || !b2 || !b3 || !b4) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	if (!m1 || !m2 || !m3 || !m4 || !m5 || !m6 || !m7) {
		fprintf(stderr, "malloc() failed\n");
		return;
	} // s(i) matrices will be used to determine the partial sums of a and b;
	int **s1, **s2, **s3, **s4, **s5, **s6, **s7, **s8, **s9, **s10;
	all_mat(&s1, sz2); all_mat(&s2, sz2); all_mat(&s3, sz2); all_mat(&s4, sz2);
	all_mat(&s5, sz2); all_mat(&s6, sz2); all_mat(&s7, sz2); all_mat(&s8, sz2);
	all_mat(&s9, sz2); all_mat(&s10, sz2);
	if (!s1 || !s2 || !s3 || !s4 || !s5 || !s6 || !s7 || !s8 || !s9 || !s10) {
		fprintf(stderr, "malloc() failed\n");
		return;
	}
	int **c1, **c2, **c3, **c4, **c5, **c6, **c7, **c8; // C = A * B
	all_mat(&c1, sz2); all_mat(&c2, sz2); all_mat(&c3, sz2); all_mat(&c4, sz2);
	all_mat(&c5, sz2); all_mat(&c6, sz2); all_mat(&c7, sz2); all_mat(&c8, sz2);
	if (!c1 || !c2 || !c3 || !c4 || !c5 || !c6 || !c7 || !c8) {
		fprintf(stderr, "malloc() failed\n");
		return;
	} // divides the matrices into portions of half the size.. block matrices
	for (int i = 0; i < sz2; ++i) {
		for (int j = 0 ; j < sz2; ++j) {
			a1[i][j] = mat1[i][j];
			a2[i][j] = mat1[i][j + sz2];
			a3[i][j] = mat1[i + sz2][j];
			a4[i][j] = mat1[i + sz2][j + sz2];
			b1[i][j] = mat2[i][j];
			b2[i][j] = mat2[i][j + sz2];
			b3[i][j] = mat2[i + sz2][j];
			b4[i][j] = mat2[i + sz2][j + sz2];
		}
	} // partial sums for calculating m1, m2..m7 from formula
	add(a1, a4, sz2, &(s1)); add(b1, b4, sz2, &(s2));
	add(a3, a4, sz2, &(s3)); subtract(b2, b4, sz2, &(s4));
	subtract(b3, b1, sz2, &(s5)); add(a1, a2, sz2, &(s6));
	subtract(a3, a1, sz2, &(s7)); add(b1, b2, sz2, &(s8));
	subtract(a2, a4, sz2, &(s9)); add(b3, b4, sz2, &(s10));
	// recursion till it gets to 1x1 matrices
	strs(s1, s2, sz2, &m1); strs(s3, b1, sz2, &m2); strs(a1, s4, sz2, &m3);
	strs(a4, s5, sz2, &m4); strs(s6, b4, sz2, &m5); strs(s7, s8, sz2, &m6);
	strs(s9, s10, sz2, &m7); // matrices are broken down into smaller ones

	add(m1, m4, sz2, &c1); subtract(c1, m5, sz2, &c2);// c3->c1(formula)
	add(c2, m7, sz2, &c3); add(m3, m5, sz2, &c4);// c4->c2
	add(m2, m4, sz2, &c5); subtract(m1, m2, sz2, &c6);// c5->c3
	add(c6, m3, sz2, &c7); add(c7, m6, sz2, &c8);// c8->c4;

	for (int i = 0; i < sz2; ++i) {
		for (int j = 0; j < sz2; ++j) {
			(*prod)[i][j] = c3[i][j]; // prod is the A * B matrix
			(*prod)[i][j + sz2] = c4[i][j]; // it's loaded up by portions
			(*prod)[i + sz2][j] = c5[i][j]; // matrix is divided in 4 blocks
			(*prod)[i + sz2][j + sz2] = c8[i][j]; // c values are from formula
		}
	} // freeing all the resources:
	free_mat(a2, sz2); free_mat(a3, sz2); free_mat(b2, sz2); free_mat(b3, sz2);
	free_mat(a1, sz2); free_mat(a4, sz2); free_mat(b1, sz2); free_mat(b4, sz2);
	free_mat(s1, sz2); free_mat(s2, sz2); free_mat(s3, sz2); free_mat(s4, sz2);
	free_mat(s5, sz2); free_mat(s6, sz2); free_mat(s7, sz2); free_mat(s8, sz2);
	free_mat(s9, sz2); free_mat(s10, sz2);
	free_mat(m1, sz2); free_mat(m2, sz2); free_mat(m3, sz2); free_mat(m4, sz2);
	free_mat(m5, sz2); free_mat(m6, sz2); free_mat(m7, sz2);
	free_mat(c1, sz2); free_mat(c2, sz2); free_mat(c6, sz2); free_mat(c7, sz2);
	free_mat(c3, sz2); free_mat(c4, sz2); free_mat(c5, sz2); free_mat(c8, sz2);
}

/* Frees the allocated memory for the container,
 * rows and columns arrays respectively 
 */
void free_memory(int ****container, int *rows, int *columns, int cnt)
{
	for (int k = 0; k < cnt; ++k) {
		for (int i = 0; i < rows[k]; ++i)
			free((*container)[k][i]);
		free((*container)[k]);
	}
	free(*container); free(rows); free(columns);
}

/* Main driver */
void main_void(char *c, int ****container, int *cnt)
{
	scanf(" %c", c);
	int nr = 1, *rows, *columns;
	allocate_memory(container, &rows, &columns, &nr);
	while (*c != 'Q') {
		if (*c == 'L') {
			if (*cnt == nr) {
				nr *= 2; // if memory is full, it doubles the size;
				reallocate_memory(container, &rows, &columns, nr);
			}
			read_matrix(container, rows, columns, cnt);
		} else if (*c == 'D') {
			command_D(*cnt, rows, columns);
		} else if (*c == 'P') {
			command_P((*container), rows, columns, *cnt);
		} else if (*c == 'C') {
			command_C(&(*container), rows, columns, *cnt);
		} else if (*c == 'M') {
			command_M(container, &rows, &columns, &(*cnt), &(nr));
		} else if (*c == 'O') {
			command_O(container, rows, columns, *cnt);
		} else if (*c == 'T') {
			command_T(container, &rows, &columns, *cnt);
		} else if (*c == 'R') {
			command_R(container, rows, columns, *cnt);
		} else if (*c == 'F') {
			int index;
			scanf("%d", &index);
			if (*cnt < index + 1 || index < 0) {
				printf("No matrix with the given index\n");
			} else {
				free_mat((*container)[index], rows[index]); //f rees matrix
				// shifts all matrix one by one, putting the freed one last;
				for (int i = index; i < *cnt - 1; ++i) {
					swap_matrices(&(*container)[i], &(*container)[i + 1]);
					swap_pos(&rows[i], &rows[i + 1]);
					swap_pos(&columns[i], &columns[i + 1]);
				}
				(*cnt)--;
				if (*cnt + 1 == nr / 2) {
					nr /= 2; // if too many were freed;
					reallocate_memory(container, &rows, &columns, nr);
				}
			}
		} else if (*c == 'S') {
			int i1, i2;
			scanf("%d %d", &i1, &i2);
			if (*cnt < i1 + 1 || *cnt < i2 + 1 || i1 < 0 || i2 < 0) {
				printf("No matrix with the given index\n");
			} else {
				if (columns[i1] != rows[i2]) {
					printf("Cannot perform matrix multiplication\n");
				} else {
					int **a;
					all_mat(&a, rows[i2]);
					strs((*container)[i1], (*container)[i2], rows[i2], &a);
					if (!a) {
						fprintf(stderr, "malloc failed\n");
						return;
					}
					if (*cnt == nr) {
						(nr) *= 2;
						reallocate_memory(container, &rows, &columns, nr);
					}
					// new matrix coords;
					rows[*cnt] = rows[i1], columns[*cnt] = columns[i2];
					// loading it up;
					all_mat(&(*container)[*cnt], rows[i1]);
					swap_matrices(&(*container)[*cnt], &a);
					free_mat(a, rows[i1]);
					(*cnt)++;
				}
			} // if command introduced is not a number:
		} else if (isdigit(*c) == 0) {
			printf("Unrecognized command\n");
		} // if any number is introduced as a command it is ignored;
		scanf(" %c", c);
	}
	free_memory(container, rows, columns, *cnt);
}

int main(void)
{
	int ***container = NULL;
	int cnt = 0;
	char command;
	main_void(&command, &container, &cnt);
	return 0;
}

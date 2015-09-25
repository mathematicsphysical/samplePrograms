#include<math.h>
#define MINEPS (1e-10)

/*
* gaussSolver 
* 部分ピボット選択付きガウスの消去法
* 
* @param[in]     n: 次元数
* @param[in]     matrix: n 行 n 列の係数行列
* @param[inout]  b: n 次元右辺ベクトル
* return: 成功 1，失敗 0
*/
int gaussSolver(int n, double *matrix, double *b){
	int i, result;
	int pivot;
	int row, col;
	double max, tmp;
	double diag, rowScale, sum;
	result = 1;

	// 前進消去を行う
	for(i = 0; i < n - 1; i++){
		// ピボット選択を行う
		pivot = i;
		max = fabs(matrix[i * n + i]);

		// 対象とする列から，最大値を探す
		for(row = i + 1; row < n; row++){
			tmp = fabs(matrix[row * n + i]);

			if(max < tmp){
				pivot = row;
				max = tmp;
			}
		}

		// 対角成分が最大値になるように，行を入れ替える
		if(pivot != i){
			for(col = 0; col < n; col++){
				tmp = matrix[i * n + col];
				matrix[i * n + col] = matrix[pivot * n + col];
				matrix[pivot * n + col] = tmp;
			}
			tmp = b[i];
			b[i] = b[pivot];
			b[pivot] = tmp;
		}
		diag = matrix[i * n + i];

		if(fabs(diag) < (double)MINEPS){
			result = 0;
			break;
		}

		for(row = i + 1; row < n; row++){
			rowScale = matrix[row * n + i] / diag;

			for(col = i; col < n; col++){
				matrix[row * n + col] -= matrix[i * n + col] * rowScale;
			}
			b[row] -= b[i] * rowScale;
		}
	}

	// 後退代入を行う
	if(result == 1){
		for(row = n - 1; row >= 0; row--){
			sum = b[row];

			for(col = row + 1; col < n; col++){
				sum -= matrix[row * n + col] * b[col];
			}
			b[row] = sum / matrix[row * n + row];
		}
	}

	return result;
}

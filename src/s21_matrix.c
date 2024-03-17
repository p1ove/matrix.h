#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = 0;
  if (rows < 1 || columns < 1) {
    error = 1;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    error = 0;
  }
  if (result->matrix != NULL) {
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = (double *)calloc(columns, sizeof(double));
      if (!result->matrix[i]) {
        for (int j = 0; j < i; j++) free(result->matrix[j]);
        free(result->matrix);
      }
    }
    error = 0;
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
  }
  A->matrix = NULL;
  A->rows = 0;
  A->columns = 0;
}

int s21_is_Emty(matrix_t *matrix) {
  int res = 0;
  if (matrix == NULL || matrix->matrix == NULL || matrix->rows <= 0 ||
      matrix->columns <= 0) {
    res = 1;
  } else {
    res = 0;
  }
  return res;
}

int s21_size_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;

  if (A->columns != B->columns) {
    res = FAILURE;
  }
  if (A->rows != B->rows) {
    res = FAILURE;
  }
  return res;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int res = SUCCESS;
  if (!(A && B)) {
    res = 0;
  } else {
    res = s21_size_eq_matrix(A, B);

    if (res == SUCCESS) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
            res = FAILURE;
          }
        }
      }
    }
  }
  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0 && s21_is_Emty(B) == 0) {
    if ((A->rows == B->rows) && (A->columns == B->columns)) {
      res = s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      res = 2;
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0 && s21_is_Emty(B) == 0) {
    if ((A->rows == B->rows) && (A->columns == B->columns)) {
      res = s21_create_matrix(A->rows, A->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      res = 2;
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0) {
    res = s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0 && s21_is_Emty(B) == 0) {
    if ((A->columns == B->rows)) {
      res = s21_create_matrix(A->rows, B->columns, result);
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < B->rows; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      res = 2;
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0) {
    res = s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[j][i];
      }
    }
  } else {
    res = 1;
  }
  return res;
}

void s21_minor(matrix_t *A, matrix_t *B, int x, int y) {
  int shift_i = 0;
  for (int i = 0; i < A->rows - 1 + shift_i; i++) {
    if (i == x) {
      shift_i++;
      continue;
    }
    int shift_j = 0;
    for (int j = 0; j < A->columns - 1 + shift_j; j++) {
      if (j == y) {
        shift_j++;
        continue;
      }
      B->matrix[i - shift_i][j - shift_j] = A->matrix[i][j];
    }
  }
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0) {
    if (A->columns == A->rows && A->columns > 1 && A->rows > 1) {
      res = s21_create_matrix(A->rows, A->columns, result);
      for (int x = 0; x < A->rows; x++) {
        matrix_t B;
        if (s21_create_matrix(A->rows - 1, A->columns - 1, &B) || res) {
          res = 1;
          s21_remove_matrix(&B);
          break;
        }
        for (int y = 0; y < A->columns && !res; y++) {
          s21_minor(A, &B, x, y);
          double det = 0;
          res = s21_determinant(&B, &det);
          if (res) {
            s21_remove_matrix(&B);
            break;
          }
          result->matrix[x][y] = det * pow(-1, x + y);
        }
        s21_remove_matrix(&B);
      }
    } else {
      res = 2;
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0) {
    if (A->columns != A->rows) {
      res = 2;
    } else {
      if (A->columns == 1 && A->rows == 1) {
        *result = A->matrix[0][0];
      } else {
        if (A->columns == 2 && A->rows == 2) {
          *result = A->matrix[0][0] * A->matrix[1][1] -
                    A->matrix[0][1] * A->matrix[1][0];
        } else {
          matrix_t B = {0};
          s21_create_matrix(A->columns - 1, A->rows - 1, &B);
          for (int y = 0; y < A->columns; y++) {
            s21_minor(A, &B, 0, y);
            double res2 = 0;
            s21_determinant(&B, &res2);
            if (!(y % 2)) {
              *result += A->matrix[0][y] * res2;
            } else {
              *result -= A->matrix[0][y] * res2;
            }
          }
          s21_remove_matrix(&B);
        }
      }
    }
  } else {
    res = 1;
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = 0;
  if (s21_is_Emty(A) == 0) {
    // res = s21_create_matrix(A->rows, A->columns, result);
    double det = 0;
    res = s21_determinant(A, &det);
    if (det && !res) {
      matrix_t B;
      s21_calc_complements(A, &B);
      matrix_t C;
      s21_transpose(&B, &C);
      s21_mult_number(&C, 1.0 / det, result);
      s21_remove_matrix(&B);
      s21_remove_matrix(&C);

    } else {
      res = 2;
    }
  } else {
    res = 1;
  }
  return res;
}

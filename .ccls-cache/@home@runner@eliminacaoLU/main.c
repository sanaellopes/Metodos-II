#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int linhaMaiorPivo(double **M, int l, int nl);
double **EliminaGauss(double **M, int *m, int *n);
void TrocaLinha(double **i, double **j);
double **CopiaMatriz(double **, int, int);
void imprimeMatriz(double **M, int l, int c);
double **alocaMatriz(int l, int c);
double **leMatriz(char *caminho, int *m, int *n);
double *alocaVetor(int n);
void imprimeSolucao(double *solucao, int n);
double *SubstiruicaoDireta(double **M, int m, double *y);
double *SubstituicaoReversa(double **M, int m, double *y);
double **multiplica(double **A, int m, int n, double **B, int p, int q);
void LU(double **, double ***, double ***, double **, int, int);

int main(int argc, char **argv) {
  double **M, **L, **U, *b, *v, ldb;
  int l, c, p, q, j;

  if (!(M = leMatriz(argv[1], &l, &c)))
    printf("\n\tNão foi possivel ler a matriz.");
  else {
    printf("\n\tMatriz lida do arquivo %s:\n", argv[1]);
    imprimeMatriz(M, l, c);
    LU(M, &L, &U, &b, l, c);
    if (L && U) {
      printf("\n\n\tMatriz L:\n");
      imprimeMatriz(L, l, l);
      printf("\n\n\tMatriz U:\n");
      imprimeMatriz(U, l, l);
      
      printf("\n\n\tSolucoes do sistema:\n");
      imprimeSolucao(SubstituicaoReversa(M,l,b), l);
    }
  }
  printf("\n");
  return 0;
}

void LU(double **M, double ***L, double ***U, double **b, int l, int c) {
  int i, j, k, p;
  double lbd, pivot;
  *L = (double**)malloc(l*sizeof(double *));
  *U = (double**)calloc(l,sizeof(double *));
  *b = (double *)malloc(l * sizeof(double));

  for(i=0;i<l;i++) (*U)[i] = (double *)malloc(l*sizeof(double));
  for(i=0;i<l;i++) (*L)[i] = (double *)calloc(l,sizeof(double));
  
  for (i = 0; i < l; i++) {
    (*b)[i] = M[i][c-1];
  }
  
  for (i = 0; i < l; i++) {
    (*L)[i][i] = 1;
    for (j = 0; j < l; j++)
      (*U)[i][j] = M[i][j];
  }

  for (j = 0; j < l - 1; j++) {
    pivot = fabs((*U)[j][j]);
    p = 0;

    for (i = j + 1; i < l; i++) {
      if (fabs((*U)[i][j]) > pivot) {
        pivot = fabs((*U)[i][j]);
        p = i;
      }
    }

    if (p) {
      TrocaLinha(&(*U)[p], &(*U)[j]);
      lbd = (*b)[p];
      (*b)[p] = (*b)[j];
      (*b)[j] = lbd;
    }

    for (i = j+1; i < l; i++) {
      lbd = (*L)[i][j] = (*U)[i][j] / (*U)[j][j];
      for (k = j+1; k < l; k++)
        (*U)[i][k] -= lbd * (*U)[j][k];
    }
  }
  printf("\n\tProduto LU:\n");
  imprimeMatriz(multiplica(*L,l,l,*U,l,l),l,l);
}

double *SubstituicaoDireta(double **L, int m, double *y) {
  int i, j;
  double soma = 0;
  double *raizes = alocaVetor(m);

  for (i = 0; i < m; i--) {
    for (j = 0; j < i - 1; j++)
      soma += L[j][i] * raizes[i];
    raizes[i] = (L[i][m - 1] - soma) / L[i][i];
    soma = 0;
  }
  return raizes;
}

double *SubstituicaoReversa(double **M, int m, double *y) {
  int i, j;
  double soma = 0;
  double *raizes = alocaVetor(m);

  for (i = m - 1; i >= 0; i--) {
    for (j = 0; j < (m - 1 - i); j++)
      soma += M[i][m - 1 - j] * raizes[m - 1 - j];
    raizes[i] = (M[i][m] - soma) / M[i][i];
    soma = 0;
  }
  return raizes;
}

double **multiplica(double **A, int m, int n, double **B, int p, int q) {
  int i, j, k;
  double **M;
  if (!(M = alocaMatriz(m, q))) {
    printf("\n\tErro de alocação em multiplica.\n");
    exit(1);
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < q; j++) {
      M[i][j] = 0;
      for (k = 0; k < n; k++)
        M[i][j] += (A[i][k] * B[k][j]);
    }
  }
  return M;
}
int linhaMaiorPivo(double **M, int l, int nl) {
  int i, linha = l;
  double maior = fabs(M[l][l]);
  for (i = l + 1; i < nl; i++) {
    if (maior < fabs(M[i][l])) {
      maior = fabs(M[i][l]);
      linha = i;
    }
  }
  return linha;
}

void imprimeSolucao(double *solucao, int n) {
  int i;
  for (i = 0; i < n; i++)
    printf("\n\t raiz[%d] = %.2lf", i + 1, solucao[i]);
}

double *alocaVetor(int n) { return (double *)malloc(n * sizeof(double)); }

double **EliminaGauss(double **M, int *m, int *n) {
  int i, j, k, k2, linha;
  double l;

  for (k = 0; k < (*n) - 1; k++) {
    linha = linhaMaiorPivo(M, k, *m);
    if (M[linha][k] == 0) {
      printf("\n\n\tSistema sem solução!\n");
      return NULL;
    }
    TrocaLinha(&M[k], &M[linha]);
    imprimeMatriz(M, *m, *n);
    for (i = k + 1; i <= (*m) - 1; i++) {
      l = M[i][k] / M[k][k];
      for (j = 0; j < *n; j++)
        M[i][j] -= (l * M[k][j]);
    }
  }

  return M;
}

void TrocaLinha(double **i, double **j) {
  double *aux;
  aux = *i;
  *i = *j;
  *j = aux;
}

double **CopiaMatriz(double **M, int l, int c) {
  double **C;
  C = alocaMatriz(l, c);
  memcpy(C, M, l * c * sizeof(double));
  return C;
}

void imprimeMatriz(double **M, int l, int c) {
  int i, j;
  printf("\n\t");
  for (i = 0; i < l; i++) {
    printf("\n\t");
    for (j = 0; j < c; j++)
      printf("%g\t", M[i][j]);
  }
}
double **alocaMatriz(int l, int c) {
  double **M;
  int i, j;
  M = (double **)malloc(l * sizeof(double *));
  for (i = 0; i < l; i++)
    M[i] = (double *)malloc(c * sizeof(double));
  return M;
}
double **leMatriz(char *caminho, int *m, int *n) {
  double **M;
  int i, j, l, c;
  FILE *arq;
  arq = fopen(caminho, "r");
  if (arq == NULL)
    return NULL;
  fscanf(arq, "%d %d", &l, &c);
  *m = l;
  *n = c;
  M = alocaMatriz(l, c);
  for (i = 0; i < l; i++)
    for (j = 0; j < c; j++)
      fscanf(arq, "%lf", &M[i][j]);
  return M;
}

#ifndef FUNCAO_H_INCLUDED
#define FUNCAO_H_INCLUDED
#include<stdbool.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define PI 3.141592653589793
#define TAM 3

bool inicializa_vetor(double **matriz);
void inicializa_matriz(double **matriz);
void MP(double ** matriz);
void MPI(double** matriz);
void MQR(double** matriz);
bool multiplica_matrizes(double **matriz, double *y, double *z);
double maior_vetor(double *z);
double* aloca_vetor();
double* proximo_Y(double *z, double alfa, int *k);
void guarda_lambda_anterior(double* lambda, double* lambda_previous);
void calcula_aproximacao(double* z, double* y, double* lambda);
bool calcula_erro(double* lambda, double *lambda_previous, double* Erro);
double retorna_modulo(double current, double previous);
bool verifica_autovalores(double* Erro, double* y, double* lambda, double P);
bool Gauss_Jordan(double** matriz, double* y, double* z);
void zera_colunas(double x[][TAM+1]);
void trocar_linhas(double x[][TAM+1], int cont, int aux);
void resolve_sistema(double A[][TAM+1], double* z);
bool verifica_autovalores_MPI(double* Erro, double* y, double* lambda, double P);
void zera_A(double** X);
void zera_A10(double** A, double U1[][TAM]);
void zera_A20(double** A, double U2[][TAM]);
void zera_A21(double** A, double U3[][TAM]);
void devolve_identidade(double x[][TAM]);
void multiplica_M(double U3[][TAM], double U2[][TAM], double U1[][TAM], double** A, double R[][TAM]);
void transposta(double U1[][TAM], double U2[][TAM], double U3[][TAM], double Q[][TAM]);
void Proximo_A(double R[][TAM], double Q[][TAM], double** A);
bool verifica_autovalores_MQR(double** A, double P);
#endif // FUNCAO_H_INCLUDED

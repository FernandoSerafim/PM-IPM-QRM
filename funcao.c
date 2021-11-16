#include "funcao.h"

bool inicializa_vetor(double **matriz){
	int i;
	//alocando os vetores linha
	for(i=0; i < TAM; i++){
		//matriz[i] = (double*)malloc(TAM*sizeof(double));
		matriz[i] = malloc(TAM*sizeof(double)); //É o mesmo efeito da linha acima
	}
return true;
}

void inicializa_matriz(double **matriz){
	int i, j;
	int auxiliary;
	auxiliary = 1;

	for(i=0; i<TAM; i++){
			for(j=0; j<TAM; j++){

				printf("Elemento da linha %d na coluna %d: ", auxiliary, j+1);
				scanf(" %lf", &matriz[i][j]);

			}
				auxiliary++;
		}
return;
}

void MP(double ** matriz){

	system("cls"); //Não é recomendado usar isso (Não me lembro do porque)

	int i;
	double P;
	double *y, *z, *lambda;
	int k = 0;
	double alfa = 0;
	double lambda_previous[TAM];
	double Erro[TAM];
	bool retorno;

	lambda = aloca_vetor();
	z = aloca_vetor();
	y = aloca_vetor();
	for(i=0;i<TAM; i++) y[i] = 1; //Y0 OU CHUTE INICIAL

	printf("Digite a precisão que deseja utilizar para o método (em decimal, por favor): ");
	scanf(" %lf", &P);

	if(multiplica_matrizes(matriz, y, z)) alfa = maior_vetor(z); //obtendo Z1 e também alfa

	y = proximo_Y(z, alfa, &k); //vamos obter Y1

	if(k>0 && multiplica_matrizes(matriz, y, z)){ //vamos obter Z2

		calcula_aproximacao(z,y, lambda); //obtenho meu lambda 1

    do
    {
        guarda_lambda_anterior(lambda, lambda_previous); //guardei lambda 1

        alfa = maior_vetor(lambda);

        y = proximo_Y(z, alfa, &k);

        multiplica_matrizes(matriz, y, z);

        calcula_aproximacao(z,y,lambda);

        if(calcula_erro(lambda, lambda_previous, Erro)) retorno = verifica_autovalores(Erro, y, lambda, P);

    }while(retorno == false);

    }


free (y);
free (z);
free (lambda);
y = NULL;		//Não é legal deixar ponteiros soltos no programa
z = NULL;		//E ajuda a evitar que hackers explorem essa falha
lambda = NULL;
}

double* aloca_vetor(){
	double *novo_vetor = (double*)malloc(TAM*sizeof(double));
return novo_vetor;
}

bool multiplica_matrizes(double **matriz, double *y, double *z){
	int i, j;
	double soma;

	for(i=0; i<TAM; i++){
		soma = 0;
		for(j=0; j<TAM; j++){
			soma = (matriz[i][j] * y[j]) + soma;
		}
		z[i] = soma;
	}
return true;
}

double maior_vetor(double *z){
	int j;
	double current, previous;
	current = 0;

	for(j=0; j<TAM; j++){

		if(z[j] < 0) previous = (-1) * z[j];
		else previous = z[j];

		if(previous > current)  current = previous;
	}

return current;
}

double* proximo_Y(double *z, double alfa, int *k){

	double* Y_K = aloca_vetor();
	int i;
	double aux;

	aux = 1.0/alfa;

	for(i=0; i<TAM;i++) Y_K[i] = aux * z[i];

*k = *k+1;
return Y_K;
}

void calcula_aproximacao(double* z, double* y, double* lambda){
	int i;
	for(i=0; i<TAM; i++)	lambda[i] = z[i]/y[i];
}

void guarda_lambda_anterior(double* lambda, double* lambda_previous){
	int i;
	for(i=0; i<TAM; i++) lambda_previous[i] = lambda[i];
return;
}

bool calcula_erro(double* lambda, double *lambda_previous, double* Erro){
	int i;
	double modulo = 0;
	double current, previous;

	for(i=0; i<TAM; i++){

		current = lambda[i];
		previous = lambda_previous[i];

		modulo = retorna_modulo(current, previous);
		Erro[i] = modulo/lambda[i];

	}
return true;
}

double retorna_modulo(double current, double previous){

	double aux;

		aux = current - previous;

		if(aux < 0){

			aux = aux * (-1);
			return aux;

		}else return aux;
}

bool verifica_autovalores(double* Erro, double* y, double* lambda, double P){

	int i, current, aux, j;
	double cont;
	double auxiliar[TAM];
	aux = 0;
	current = 0;

	for(i=0; i<TAM; i++){

		if(Erro[i] > P)	aux++;

	}

	if(aux == TAM)	return false;

	for(i=0; i<TAM; i++){

	if(Erro[i] < P){

		auxiliar[i] = lambda[i];

		}else{
			auxiliar[i] = 10000.0;
		}

	}

	for(j=0; j<TAM; j++){

		for(i=0; i<TAM; i++){

			if(auxiliar[i] < auxiliar[j]){

				cont = auxiliar[i];
				auxiliar[i] = auxiliar[j];
				auxiliar[j] = cont;

			}
		}
	}



	for(i=0; i<TAM; i++){
		if(auxiliar[i] != 10000.0)	printf("\nO NOSSO %d AUTOVALOR É: %5f", current, auxiliar[i]);
		current++;
	}
			printf("\n");
			printf("\nCOM RESPECTIVO AUTOVETOR:\n\n");

			for(j=0; j<TAM; j++)	printf("%5f\n", y[j]);

return true;
}

void MPI(double** matriz){

	system("cls");

	int i;
	double P;
	double *y, *z, *lambda;
	int k = 0;
	double alfa = 0;
	double lambda_previous[TAM];
	double Erro[TAM];
	bool retorno;

	lambda = aloca_vetor();
	z = aloca_vetor();
	y = aloca_vetor();
	for(i=0;i<TAM; i++) y[i] = 1; //Y0

	printf("Digite a precisão que deseja utilizar para o método (em decimal, por favor): ");
	scanf(" %lf", &P);

	if(Gauss_Jordan(matriz, y, z)) alfa = maior_vetor(z); //obtendo Z1 e também alfa

	y = proximo_Y(z, alfa, &k); //vamos obter Y1

	if(k>0 && Gauss_Jordan(matriz, y, z)){ //vamos obter Z2

		calcula_aproximacao(z,y, lambda); //obtenho meu lambda 1

			do
			{
					guarda_lambda_anterior(lambda, lambda_previous); //guardei lambda 1

					alfa = maior_vetor(z);

					y = proximo_Y(z, alfa, &k);

					Gauss_Jordan(matriz, y, z);

					calcula_aproximacao(z,y,lambda);

				if(calcula_erro(lambda, lambda_previous, Erro)) retorno = verifica_autovalores_MPI(Erro, y, lambda, P);

			}while(retorno == false);

}

free (y);
free (z);
free (lambda);
y = NULL;
z = NULL;
lambda = NULL;
}

bool Gauss_Jordan(double** matriz, double* y, double* z){

/*A MATRIZ A É O SISTEMA DE EQUAÇÕES QUE POSSUI A PARTE QUADRADA IGUAL A NOSSA
MATRIZ A QUE É LIDA PELO TECLADO E A ÚLTIMA COLUNA IGUAL AO NOSSO VETOR Y, QUE
INICIALMENTE, TEM COMO CHUTE INICIAL [1,1,1] PARA O SISTEMA 3X3*/

	double A[TAM][TAM+1];
	int i, j;

    /*MINHA MATRIZ "A" RECEBE OS COEFICIENTES DA MATRIZ PRINCIPAL. FIZ ISSO PORQUE
    NESSA ITERAÇÃO NÃO DESEJO ALTERAR OS VALORES DELA. PORTANTO, FIZ UMA CÓPIA*/
    for(i=0; i<TAM; i++){
		for(j=0; j<TAM; j++){
                A[i][j] = matriz[i][j];
		}
	}

    for(i=0; i<TAM; i++){
            A[i][TAM] = y[i];/*A ÚLTIMA COLUNA DE A RECEBE PS TERMOS INDEPENDENTES*/
	}

    zera_colunas(A);
    resolve_sistema(A, z);
	return true;

}

void zera_colunas(double x[][TAM+1]){

    int i, j, sistema_SPI = 0;
    int cont = 0, aux = 0;
    double constante = 0.0;

  do
  {
      if(x[cont][cont] == 0 && x[aux+1][cont] != 0) trocar_linhas(x, cont, aux);

    for(i=cont; i<TAM; i++)
        {
          constante =  (x[aux + 1][cont] / x[cont][cont]) * -1;

            for(j=0; j<TAM+1; j++)
                {
                    if(i == cont) break;
                    else x[i][j] = (x[cont][j] * constante) + x[i][j];
                }
            if(i != cont ) aux++;
        }

    cont++;
    aux = cont;
  }while(cont != TAM-1);

 /*FOR QUE VERIFICA SE TODOS OS ELEMENTOS DA ÚLTIMA LINHA SÃO NULOS*/
    for(i=0; i<TAM-1; i++)
    {
        if(x[TAM-1][i] == 0) sistema_SPI++;
    }

    if(sistema_SPI == TAM)
    {
        printf("O sistema possui infinitas soluções\n");
        exit(0);
    }

}

void trocar_linhas(double x[][TAM+1], int cont, int aux){
    int j;
    double recebe_valor;

    //CONT = previous line       ---->             x[cont][cont]
    //AUX  = current line        ---->             x[aux+1][cont]

    for(j=0; j<TAM+1; j++)
    {
        recebe_valor = x[aux + 1][j];
        x[aux + 1][j] = x[cont][j];
        x[cont][j] = recebe_valor;
    }
}

void resolve_sistema(double x[][TAM+1], double* z){

    //double respostas[TAMANHO];
    double soma = 0.0;
    int i, j, aux = 0;
    int elementos = 1, contador = 0;
    int linha = TAM - 2;
    int coluna = TAM - 1; //NÚMERO COLUNAS COM ELEMENTOS A SEREM SOMADOS

    for(i = 0; i < TAM; i++) z[i] = PI;

    z[TAM-1] = x[TAM-1][TAM]/x[TAM-1][TAM-1];

    /*FOR PARA MULTIPLICAR ELEMENTOS QUE VÃO TROCAR DE SINAL*/
    for(i = 0; i < TAM ; i++)
    {
        for(j = elementos; j < TAM ; j++)
            {
                x[i][j] = x[i][j] * -1;
            }
        elementos++;
    }

do
{
    /*FOR PARA DETERMINAR QUAL POSIÇÃO DA MATRIZ JÁ POSSUI RESPOSTA*/
    for(i = 0; i < TAM ; i++)
    {
            if(z[i] != PI)
            {
                aux = i;
                break;
            }
    }

   /*LOOPING QUE MULTIPLICA TODA A COLUNA QUE TIVER COM SUA VARIÁVEL JÁ ENCONTRADA*/
    do
    {
        x[contador][aux] = x[contador][aux] * z[aux];
        contador++;

    }while(contador != TAM);

    /*FOR PARA REALIZAR O SOMATÓRIO DA LINHA E OBTER A NOVA RESPOSTA*/
    for(j = coluna; j < TAM + 1; j++)
        {
                soma =   x[linha][j] + soma;
                if(j == TAM) soma = soma/x[linha][linha];
        }
z[linha] = soma;
coluna = coluna - 1;
linha = linha - 1;
soma = 0.0;
contador = 0;
}while(linha >= 0);


}

bool verifica_autovalores_MPI(double* Erro, double* y, double* lambda, double P){
	int i, aux, j;
	double cont, current;
	double auxiliar[TAM];
	aux = 0;


	for(i=0; i<TAM; i++){

		if(Erro[i] > P)	aux++;

	}

	if(aux == TAM) return false;


	for(i=0; i<TAM; i++){

	if(Erro[i] < P){

		auxiliar[i] = lambda[i];

		}else{
			auxiliar[i] = 100000.0;
		}

	}



	for(j=0; j<TAM; j++){

		for(i=0; i<TAM; i++){

			if(auxiliar[i] < auxiliar[j]){

				cont = auxiliar[i];
				auxiliar[i] = auxiliar[j];
				auxiliar[j] = cont;

			}
		}
	}

	for(i=0; i<TAM; i++){
		if(auxiliar[i] != 100000.0)	current = auxiliar[i];
		break;
	}

			printf("\nO NOSSO AUTOVALOR É: %5f", 1.0/current);
			printf("\n");
			printf("\nCOM RESPECTIVO AUTOVETOR:\n\n");
			for(j=0; j<TAM; j++){
				printf("%5f\n", y[j]);
			}

return true;
}

void MQR(double** matriz){

	system("cls");

	double U1[TAM][TAM], U2[TAM][TAM], U3[TAM][TAM], R[TAM][TAM];
	double Q[TAM][TAM];
	int i, j;
	double precisao;
    double **A;
	bool retorno;

	A = (double**)malloc(TAM * sizeof(sizeof(double)));

	if(inicializa_vetor(A)) zera_A(A);

	for(i=0; i<TAM; i++){
		for(j=0; j<TAM; j++) A[i][j] = matriz[i][j];
	}

	printf("Digite a precisão que deseja utilizar para o método (em decimal, por favor): ");
	scanf(" %lf", &precisao);

do{
	devolve_identidade(U1);
	if(A[1][0] != 0)	zera_A10(A, U1); //obtém matriz para zerar a10 (se preciso for)

	devolve_identidade(U2);
	if(A[2][0] != 0)	zera_A20(A, U2); //obtém matriz para zerar a20 (se preciso for)

	devolve_identidade(U3);
	if(A[2][1] != 0)	zera_A21(A, U3); //obtém matriz para zerar a21 (se preciso for)

	multiplica_M(U3, U2, U1, A, R); //obtém R

	transposta(U1, U2, U3, Q); //obtém Q

	Proximo_A(R, Q, A); //obtém próximo A

	retorno = verifica_autovalores_MQR(A, precisao); //verifica os autovalores

}while(retorno == false);

	printf("\n\nMatriz 'A' FINAL\n");
	for(i=0; i<TAM; i++){
	for(j=0; j<TAM; j++) printf("%5f ", A[i][j]);
	printf("\n");
	}

free(A);
}

void zera_A(double** A){
	int i, j;

	for(i=0; i<TAM; i++){
		for(j=0; j<TAM; j++)	A[i][j] = 0.0;
	}
return;
}

void zera_A10(double** A, double U1[][TAM]){

	int i, j;
	int q, p;

	q = 1;
	p = 0;

	//LÓGICA PARA O COSSENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( (i == q && j == q) || (i == p && j == p) ){

				U1[i][j] = A[p][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

		}

	}

	//LÓGICA PARA O SENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( i == p && j == q ){

				U1[i][j] = A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

			if( i == q && j == p ){

				U1[i][j] = -1 * (A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) ) );

			}
		}
	}
}

void zera_A20(double** A, double U2[][TAM]){

	int i, j;
	int q, p;

	q = 2;
	p = 0;

	//LÓGICA PARA O COSSENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( (i == q && j == q) || (i == p && j == p) ){

				U2[i][j] = A[p][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

		}

	}

	//LÓGICA PARA O SENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( i == p && j == q ){

				U2[i][j] = A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

			if( i == q && j == p ){

				U2[i][j] = -1 * (A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) ) );

			}

		}

	}

}

void zera_A21(double** A, double U3[][TAM]){

	int i, j;
	int q, p;

	q = 2;
	p = 1;

	//LÓGICA PARA O COSSENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( (i == q && j == q) || (i == p && j == p) ){

				U3[i][j] = A[p][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

		}

	}

	//LÓGICA PARA O SENO
	for(i = 0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			if( i == p && j == q ){

				U3[i][j] = A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) );

			}

			if( i == q && j == p ){

				U3[i][j] = -1 * (A[q][p]/( sqrt( pow(A[p][p],2) + pow(A[q][p],2) ) ) );

			}

		}

	}

}

void devolve_identidade(double x[][TAM]){
int i, j;

		for(i=0; i<TAM; i++){

			for(j=0; j<TAM; j++){

				if(j != i) x[i][j] = 0;
				else x[i][j] = 1;

			}

		}
}

void multiplica_M(double U3[][TAM], double U2[][TAM], double U1[][TAM], double** A, double R[][TAM]){

	int i, j, cont;
	float aux[TAM][TAM], guarda[TAM][TAM], soma;

	cont = 0;
	soma = 0.0;

for(i=0; i<TAM; i++){
		for(j=0; j<TAM; j++){
			aux[i][j] = 0.0;
			guarda[i][j] = 0.0;
		}
	}
//MULTIPLICA U3 E U2 E GUARDA EM AUX
do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  (U3[cont][j] *	U2[j][i]) + soma;

		}

		aux[cont][i] = soma;

	}

	cont++;

}while(cont != 3);

cont = 0;
//MULTIPLICA aux e U1 E GUARDA EM guarda
do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  (aux[cont][j] *	U1[j][i]) + soma;

		}

		guarda[cont][i] = soma;

	}

	cont++;

}while(cont != 3);

cont = 0;
//MULTIPLICA guarda e A E GUARDA EM R
do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  (guarda[cont][j] *	A[j][i]) + soma;

		}

		R[cont][i] = soma;

	}

	cont++;

}while(cont != 3);


}

void transposta(double U1[][TAM], double U2[][TAM], double U3[][TAM], double Q[][TAM]){

	int i, j;
	float vetor1[TAM][TAM], vetor2[TAM][TAM], vetor3[TAM][TAM];
	float armazena[TAM][TAM];
	float soma;
	int cont;
	cont = 0;

//FAZENDO A TRANSPOSTA DE CADA "U"
	for(i=0; i<TAM; i++){

		for(j=0; j<TAM; j++){

			vetor1[i][j] = U1[j][i];
			vetor2[i][j] = U2[j][i];
			vetor3[i][j] = U3[j][i];
			armazena[i][j] = 0.0;

		}
	}
//AGORA QUE JÁ TENHO AS TRANSPOSTAS, POSSO MULTIPLICÁ-LAS
do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  (vetor1[cont][j] *	vetor2[j][i]) + soma;

		}

		armazena[cont][i] = soma;

	}

	cont++;

}while(cont != 3);

cont = 0;
do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  (armazena[cont][j] * vetor3[j][i]) + soma;

		}

		Q[cont][i] = soma;

	}

	cont++;

}while(cont != 3);

}

void Proximo_A(double R[][TAM], double Q[][TAM], double** A){

	int i, j, cont;
	float soma;

	cont = 0;

do{

	for(i=0; i<TAM; i++){

		soma = 0.0;

		for(j=0; j<TAM; j++){

			soma =  ( R[cont][j] * Q[j][i] ) + soma;

		}

		A[cont][i] = soma;

	}

	cont++;

}while(cont != 3);
}

bool verifica_autovalores_MQR(double** A, double P){
    int i, j;
	if(A[1][0] < P && A[2][0] < P && A[2][1] < P){
		for(i=0; i<TAM; i++){
			for(j=0; j<TAM; j++){
				if(j == i) printf("\nAUTOVALOR: %5f\n\n", A[i][j]);
			}
		}
		return true;
	}else return false;
}

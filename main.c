#include "funcao.h"
#include "funcao.c"
#include<locale.h>

int main(int argc, char **argv){

	setlocale(LC_ALL, "Portuguese_Brazil");
	printf("Favor usar v�rgula ao inv�s de ponto quando for designar n�meros decimais!\n\n\n");

	double **matriz;
	int opt, i;

	//matriz = (double**)malloc(TAM * sizeof(sizeof(double)));
	matriz = malloc(TAM * sizeof(double*));	 //O efeito � igualzinho ao de cima!
	if(inicializa_vetor(matriz)) inicializa_matriz(matriz);

		do{
				printf("\n\nDigite o m�todo que deseja calcular autovalores\n");
				printf("[1] M�todo das pot�ncias\n[2]M�todo das pot�ncias inversas\n[3]M�todo QR\n");
				scanf(" %d", &opt);

			switch(opt){
				case 0:
					printf("ENCERRANDO SISTEMA\n");
					break;
				case 1:
					MP(matriz);
					break;
				case 2:
					MPI(matriz);
					break;
				case 3:
					MQR(matriz);
					break;
				default:
					printf("Op��o inv�lida\n");
				}

	}while(opt != 0);

for(i=0; i<TAM; i++){
	free (matriz[i]);
	matriz[i] = NULL;
}
free(matriz); matriz = NULL;
return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TRUE 1
#define FALSE 0
typedef struct _cidade_{
    int x;
    int y;
    int dema;
}city;
city *citys;
typedef struct _sol_{
    float dist;
    int capacidade;
    int demanda;
    int car;
    int tam;
    int *clientes;
}solucao;
solucao *sol;
float **custo;
unsigned int chamada=0;
void free_tudo(solucao *s,int tam){
    int i;
    for(i =0 ;i < tam ;i++){
        free(s[i].clientes);
    }
    free(s);
}
void free_matriz(float **m,int tam){
    int i;
    for(i = 0 ; i < tam ;i++){
        free(m[i]);
    }
    free(m);
}
float **aloca_matriz(int tam){
    float **m;
    int i;
    m=(float **)malloc(tam*sizeof(float*));
    if(m!=NULL){
        for(i=0;i<tam;i++){
            m[i] = (float *)malloc(tam*sizeof(float));
            if(m[i] == NULL){
                return NULL;
            }
        }
        return m;
    }
    return NULL;
}
solucao *aloca_sol(int tam,int tam_client,int cap){
    int i;
    solucao *s;
    s = (solucao *) malloc(tam*sizeof(solucao));
    if(s != NULL){
        for(i = 0 ; i<tam;i++){
            s[i].clientes = (int *)malloc(tam_client*sizeof(int));
            s[i].capacidade = cap;
            s[i].car=0;
            s[i].demanda=0;
            s[i].dist=0.0;
            s[i].tam = 0;
            if(s[i].clientes == NULL){
                return NULL;
            }
        }
        return s;
    }
    return NULL;
}
int calc_dema(int v[],int tam){
    int i,dem = 0;
    if(sizeof(v)<=0){
        printf("Nao factivel,sugestao diminua a capacidade,solucao saturada.%d",sizeof(v));
        return;
    }
    for(i=0 ; i< tam ;i++){
        dem+=citys[v[i]].dema;
    }
    return dem;
}
void copia_solu(solucao *s,solucao *s2,int car){
    int i,j;
    if(s == NULL) return;
    for(i = 0 ; i < car ; i ++){
        s2[i].capacidade = s[i].capacidade;
        s2[i].car = s[i].car;
        s2[i].demanda = s[i].demanda;
        s2[i].dist = s[i].dist;
        s2[i].tam = s[i].tam;
        for(j=0;j<s[i].tam;j++){
            s2[i].clientes[j] = s[i].clientes[j];
        }
    }
}
int isFULL(int v[],int tam){
    int i;
    for(i=0 ; i < tam ; i++){
        if(v[i] != -1){
            return FALSE;
        }
    }
    return TRUE;
}
void printaVet(int v[],int tam){
    int i;
    for(i = 0;i<tam;i++){
        printf("%d ",v[i]);
    }
    printf("\n");
}
void zerar(int v[],int tam){
int i = 0;
for(i = 0 ; i < 5 ; i ++) v[i] = -1;
}
int isIn(int v[],int tam,int k){
    int i;
    for(i = 0 ; i < tam ;i++){
        if(v[i] == k){
            return i;
        }
    }
    return -1;
}
int maior(float v[],int tam,int visit[]){
    float m = -10000000.0;
    int i,j=0;
    for( i = 0 ; i < tam ; i++){
        if(visit[i]!=-1){
            if(m < v[i]){
                m = v[i];
                j = i;
            }
        }
    }
    return j;
}
int menor(float v[],int tam,int visit[]){
    float m = 10000000000000000000000.00;
    int i,j=0;
    for( i = 0 ; i < tam ; i++){
        if(v[i]!=0.0 && visit[i]!=-1){
            if(m > v[i]){
                m = v[i];
                j = i;
            }
        }
    }
    return j;
}
float calcula_distancia(int xa,int ya,int xb,int yb){
    double a,b;
    float r;
    a = xa-xb;
    b = ya-yb;
    r = (float)sqrt(pow(a,2)+pow(b,2));
    return r;
}
int ler_Arquivo(FILE *a,FILE *b){
    int aux,aux2,aux3;
    do{
        fscanf(a,"%d %d",&aux,&aux2);
        if(!feof(a)){
            citys = (city *)realloc(citys,aux*sizeof(city));
            citys[aux-1].dema = aux2;//no 1 == 0
        }
    }while(!feof(a));
    do{
        fscanf(b,"%d %d %d",&aux,&aux2,&aux3);
        citys[aux-1].x = aux2;
        citys[aux-1].y = aux3;
    }while(!feof(b));
    fclose(a);
    fclose(b);
    return aux;
}
float calc_custo(solucao *s,int car){
    float c=0.0;
    int i;
    chamada++;
    if(s!=NULL){
        for(i = 0 ; i < car;i++){
            c+=s[i].dist;
        }
    }else{
        printf("Erro");
    }
    return c;
}
int efactivel(solucao *s,int car,int tam){
    int i,j,k,l,flag;
    for(i = 0 ; i < car ; i++){
        if(s[i].demanda > s[i].capacidade){
            return FALSE;
        }
    }
    for( i = 0 ; i < tam ;i++){
        flag = FALSE;
        for( j = 0 ; j< car ;j++){
            for(k = 0 ;k < s[j].tam ;k++){
                if(s[j].clientes[k] == i){
                    flag = TRUE;
                }
            }
        }
        if(flag == FALSE){
            return FALSE;
        }
    }
    for( i = 0 ; i < car ; i++){
        for(j = 1 ; j < s[i].tam ; j++){
            for(k = j+1;k < s[i].tam;k++){
                if(s[i].clientes[j] == s[i].clientes[k]){
                    return FALSE;
                }
            }
        }
        for(j = i+1 ; j < car ; j++){
            for( k = 1 ; k < s[i].tam ; k++){
                for(l = 1 ; l < s[j].tam ; l++){
                    if(s[i].clientes[k] == s[j].clientes[l]){
                        return FALSE;
                    }
                }
            }
        }
    }
    return TRUE;
}
solucao * constroiSolucaoInicialGulosa(int car,int cap,int t,float alpha){
    int i,j,k,l,n,tam=t,nCand,jav=0,aux,flag;
    int jaVisitado[tam],lrc[tam],demanda=0;
    float teto;
    float me,ma;
    sol = aloca_sol(car,tam,cap);
    custo = aloca_matriz(tam);
    if(custo==NULL && sol == NULL){
        printf("erro ao alocar matriz");
        free_matriz(custo,tam);
        return;
    }
    for(j = 0;j < tam;j++){
       for(k = 0;k<tam;k++){
            custo[j][k]=calcula_distancia(citys[j].x,citys[j].y,citys[k].x,citys[k].y);
        }
    }
    for(i =0;i<tam;i++){
        jaVisitado[i]=0;
    }
    jaVisitado[0] = -1;
    jav=0;
    k=0;
    while(isFULL(jaVisitado,tam) == FALSE && jav < car){
        demanda=0;
        i=0;
        k=0;
        sol[jav].clientes[0]=0;
        k++;
        while(demanda < cap){
            nCand = 0;
            me=custo[i][menor(custo[i],tam,jaVisitado)];
            ma=custo[i][maior(custo[i],tam,jaVisitado)];
            teto = me+alpha*(ma-me);
            for(j=0;j<tam;j++){
                if(jaVisitado[j]!=-1 ){
                    if(custo[i][j] <=teto){
                        lrc[nCand] = j;
                        nCand++;
                    }
                }
            }
            if(nCand == 0){
                break;}
            if((demanda+citys[i].dema)<cap){
                i=lrc[random()%nCand];
                jaVisitado[i]=-1;
                demanda+=citys[i].dema;
                sol[jav].clientes[k] = i;
                k++;
            }else{
                break;
            }
        }
        sol[jav].capacidade = cap;
        sol[jav].car = jav;
        sol[jav].demanda = demanda;
        sol[jav].tam = k;
        jav++;
    }
    if(isFULL(jaVisitado,tam) == FALSE){
        for(i = 1 ; i < tam ;i++){
            for(j = 0 ; j < car ; j++){
                if(jaVisitado[i] != -1){
                    sol[j].clientes[sol[j].tam] = i;
                    sol[j].demanda+=citys[i].dema;
                    sol[j].tam++;
                    if((sol[j].demanda+calc_dema(sol[j].clientes,sol[j].tam))<=cap){
                        sol[j].demanda=calc_dema(sol[j].clientes,sol[j].tam);
                        jaVisitado[i] = -1;
                    }else{
                        sol[j].tam--;
                        sol[j].demanda-=citys[i].dema;
                    }
                }
            }
        }
        for(j=0;j<car;j++){
            if(sol[j].demanda < sol[j].capacidade){
                for( i = 0 ; i < tam ;i++){
                    if(jaVisitado[i]!=-1){
                        sol[j].clientes[sol[j].tam]=i;
                        sol[j].demanda+=citys[i].dema;
                        sol[j].tam++;
                        jaVisitado[i] = -1;
                    }
                }
            }
        }
        for(i = 0 ; i < car ; i++){
            for(j = 0 ; j < car ;j++){
                for( k = 1 ;k<sol[i].tam;k++){
                    for(l = 1;l<sol[j].tam;l++){
                       aux = sol[i].clientes[k];
                       sol[i].clientes[k] = sol[j].clientes[l];
                       sol[j].clientes[l] = aux;
                       if((calc_dema(sol[i].clientes,sol[i].tam)<= cap)&&(calc_dema(sol[j].clientes,sol[j].tam)<= cap)){
                            sol[i].demanda = calc_dema(sol[i].clientes,sol[i].tam);
                            sol[j].demanda = calc_dema(sol[j].clientes,sol[j].tam);
                       }else{
                            aux = sol[i].clientes[k];
                            sol[i].clientes[k] = sol[j].clientes[l];
                            sol[j].clientes[l] = aux;
                       }
                    }
                }
            }
        }
    }
    for( i = 0 ; i < car ; i++){
        sol[i].dist=0.0;
        for( j = 0 ; j < sol[i].tam-1 ; j++){
            sol[i].dist+=custo[sol[i].clientes[j]][sol[i].clientes[j+1]];
        }
        //sol[i].dist+=custo[sol[i].clientes[sol[i].tam-1]][0];
    }
    return sol;
}
solucao* gerar_vizinho(solucao* solu,int car,int tam,int cap){
    int i = random()%car;
    int j;
    int k,l,n,m,aux;
    solucao *s = aloca_sol(car,tam,cap);
    if(s != NULL && solu[i].tam>0 && solu!=NULL){
        j = random()%solu[i].tam;
        for(k = 0 ; k < car ; k++){
            s[k].capacidade = solu[k].capacidade;
            s[k].car = solu[k].car;
            s[k].demanda = solu[k].demanda;
            s[k].dist = solu[k].dist;
            s[k].tam = solu[k].tam;
            for(l = 0 ; l < s[k].tam ;l++){
                s[k].clientes[l] = solu[k].clientes[l];
            }
        }
        if(j == 0){ j++;}
        for(k = 0 ; k < car ; k++){
            if(i!=k){
                for(n = 1 ; n < s[k].tam;n++){
                    aux = s[i].clientes[j];
                    s[i].clientes[j] = s[k].clientes[n];
                    s[k].clientes[n] = aux;
                    if(calc_dema(s[k].clientes,s[k].tam)<cap){
                        s[k].demanda = calc_dema(s[k].clientes,s[k].tam);
                        s[i].demanda = calc_dema(s[i].clientes,s[i].tam);
                    }else{
                        aux = s[i].clientes[j];
                        s[i].clientes[j] = s[k].clientes[n];
                        s[k].clientes[n] = aux;
                        s[k].demanda = calc_dema(s[k].clientes,s[k].tam);
                        s[i].demanda = calc_dema(s[i].clientes,s[i].tam);
                    }
                }
            }
        }
        for(i = 0 ; i < car ; i++){
            s[i].dist=0.0;
            for(j=0;j<s[i].tam-1 ; j++){
                s[i].dist += custo[s[i].clientes[j]][s[i].clientes[j+1]];
            }
            //s[i].dist+=custo[s[i].clientes[s[i].tam-1]][0];
        }
    }else{
        return NULL;
    }
    return s;
}
solucao* SA(int car,int tam,float alpha,float temperatura_ini,float temperartura_final,int maxInter,int cap,solucao *s0){
    solucao *vizinho=NULL;
    solucao *best=NULL,*t=NULL;
    float temp = temperatura_ini;
    int inter=0;
    float delta=0.0;
    float aceita=0.0;
    best = aloca_sol(car,tam,cap);
    vizinho = aloca_sol(car,tam,cap);
    if(best!=NULL && s0!=NULL){
        copia_solu(s0,best,car);
    }else{
        printf("erro null");
        return NULL;
    }
    int i,j;
    while(temp>temperartura_final){
        inter=0;
        while(inter < maxInter){
            while(t == NULL){
                t=gerar_vizinho(best,car,tam,cap);
                copia_solu(t,vizinho,car);
            }
            delta = calc_custo(vizinho,car) - calc_custo(best,car);

            if(delta <=0.0){
               copia_solu(vizinho,best,car);
            }else{
                aceita =((float)(random())/(float)(RAND_MAX));
                if(aceita < exp(-delta/temp)){
                    copia_solu(vizinho,best,car);
                }
            }
            inter++;
        }
        temp = alpha*temp;
    }
    return best;
}
void grasp(int Int_max,int car,int tam,int cap){
    solucao *s1,*s2,*f=NULL;
    int i = 0,j;
    float temp;
    while(i < Int_max){
        do{
            //0.2 alpha cmin>=c<=cmin+alpha*(max-min)
            s2=constroiSolucaoInicialGulosa(car,cap,tam,0.2);
            //alpha taxa de resfriamento a¢[0,1]
            s1 = SA(car,tam,0.9,10,0.002,100,cap,s2);
        }while(efactivel(s1,car,tam)==FALSE && efactivel(s2,car,tam)==FALSE);
        if(s1 !=NULL && s2 !=NULL){
            if(f==NULL){
                if(calc_custo(s1,car) < calc_custo(s2,car)){
                    f = s1;
                }else{
                    f = s2;
                }
            }else{
                if(calc_custo(s1,car) <= calc_custo(s2,car)){
                    if(calc_custo(f,car)<=calc_custo(s1,car)){
                        f = s1;
                    }
                }else{
                    if(calc_custo(f,car)<=calc_custo(s2,car)){
                        f = s1;
                    }
                }
            }
        }
        i++;
    }
    for(i = 0 ; i < car ; i++){
        printf("%d :carro\n",i);
        for( j = 0 ; j < f[i].tam ; j++){
            printf("%d ",f[i].clientes[j]);
        }
        printf("0 %d <-dema\n",f[i].demanda);
    }
    printf("%6.3f <-custo %u <-numero de chamadas a funcao objetivo\n",calc_custo(f,car),chamada);
}
int main()
{
    citys = (city*)malloc(sizeof(city*));
    FILE *arq = fopen("demanda.vrp","r");
    FILE *arq_d = fopen("nos.vrp","r");
    FILE *c = fopen("carro.vrp","r");
    int tam,carros,capacidade,i,j;
    float temp;
    solucao *s;
    if(citys != NULL && arq!=NULL && arq_d!=NULL && c!=NULL){
        srand((unsigned)time(NULL));
        tam = ler_Arquivo(arq,arq_d);
        fscanf(c,"%d %d",&carros,&capacidade);
        printf("clientes:%d carros:%d capacidade:%d\n",tam,carros,capacidade);
        sol = aloca_sol(carros,tam,capacidade);
        grasp(10,carros,tam,capacidade);
    }
    free_matriz(custo,tam);
    free_tudo(sol,carros);
    free(citys);
    return 0;
}

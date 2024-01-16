#include<iostream>
#include<mpi.h>

int sumar(int* tmp, int n){
    int suma=0;
    for(int i =0;i<n;i++){
        suma=suma+tmp[i];
    }
    return suma;
}


int main(int argc, char** argv){

    MPI_Init(&argc,&argv);

    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
  //  std::printf("Rank %d of %d process \n",rank,nprocs);
    
    int data[100];
    
    if(rank == 0){
        std::printf("Total de ranks: %d", nprocs);

        for(int i=0;i<100;i++){
            data[i]=i;
        }


        for(int rankid=1;rankid<nprocs; rankid++){
            std::printf("\n Rank_0 Enviando datos a rank_%d \n",rankid);
            int start= (rankid*25);
            MPI_Send(&data[start],25,MPI_INT,rankid,0,MPI_COMM_WORLD);
        }

        int suma_parcial = sumar(data,25);
        
        int suma_ranks[4];
        suma_ranks[0] = sumar(data,25);
        for(int rankid=1;rankid<nprocs; rankid++){
            MPI_Recv(&suma_ranks[rankid],1,MPI_INT,rankid,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        }

        std::printf("sumas parciales: %d, %d, %d, %d\n", suma_ranks[0],suma_ranks[1],suma_ranks[2],suma_ranks[3]);
        int suma_total=sumar(suma_ranks,4);
        
        std::printf("Suma total = %d \n", suma_total);

    }else 
    {
        std::printf(" Rank %d Recibiendo datos \n",rank);
        
        MPI_Recv(data,25, MPI_INT, 0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        
        
        //std::string str="";
        //for(int i=0;i<10;i++){
        //    str = str + std::to_string(data[i])+",";
        //}
        //std::printf("\n Rank_%d recibindo datos =====> %s \n",rank,str.c_str());
    
        int suma_parcial = sumar(data,25);

        std::printf("Rank_%d enviando suma parcial %d\n", rank, suma_parcial);

        MPI_Send(&suma_parcial, 1, MPI_INT,0,0,MPI_COMM_WORLD);    
    
    }




    MPI_Finalize();

    return 0;
}
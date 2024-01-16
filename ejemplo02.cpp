#include<iostream>
#include<mpi.h>

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


        for(int i=1;i<nprocs; i++){
            std::printf("\n Rank_0 Enviando datos a rank_%d \n",i);
            MPI_Send(
                data, //datos
                100, //cuantos datos
                MPI_INT, //tipo de dato
                i, // rank destino
                0, //tag
                MPI_COMM_WORLD); //grupo

        }


    }else if(rank==1){

        std::printf("\nRank %d Recibiendo datos ",rank);

        MPI_Recv(
            data, //datos
            100, //cuantos
            MPI_INT, //tipo de dato
            0, //rank - origen
            0, //tag
            MPI_COMM_WORLD, //grupo
            MPI_STATUS_IGNORE); //status

        std::string str="";
        
        for(int i=0;i<10;i++){
            //std::printf("%d, ",data[i]);
            str = str + std::to_string(data[i])+",";
                    std::printf("\nRank_%d recibiendo datos ===> %d",rank,data[i]);
        }

    }




    MPI_Finalize();

    return 0;
}
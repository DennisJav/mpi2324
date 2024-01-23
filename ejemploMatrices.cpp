#include<iostream>
#include<mpi.h>
#include<cmath>
#include<vector>

#define MATRIX_DIMENSION 25

void matrix_multi(double * A, double * B, double * C, int rows, int cols){
    for(int i=0;i<rows;i++){
        double tmp=0;
        for(int j=0; j<cols;j++){
            tmp=tmp+A[i*cols+j];
        }
    }
}


int main(int argc, char** argv){

    MPI_Init(&argc,&argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int rows_per_rank;
    int rows_alloc = MATRIX_DIMENSION;
    int padding;

    if(MATRIX_DIMENSION%nprocs!=0){
            rows_alloc=std::ceil((double)MATRIX_DIMENSION/nprocs)*nprocs;
            padding = rows_alloc - MATRIX_DIMENSION;
    }

    rows_per_rank = rows_alloc / nprocs;

    if(rank==0){
        //----Imprimir  informacion 
        std::printf("Dimension: %d, rows_alloc: %d, rows_per_rank: %d, padding: %d \n",
                    MATRIX_DIMENSION,rows_alloc,rows_per_rank,padding);
        
        //------
        std::vector<double> A(MATRIX_DIMENSION*rows_alloc);
        std::vector<double> B(MATRIX_DIMENSION);
        std::vector<double> C(rows_alloc);

        for(int i =0; i<MATRIX_DIMENSION*MATRIX_DIMENSION; i++){
            for(int j=0;j<MATRIX_DIMENSION;j++){
                int index=i*MATRIX_DIMENSION+j;
                A[index]=i;
            }
        
        }


        for(int i =0; i<MATRIX_DIMENSION; i++){
            B[i]=1;
        }

        //--------Enviar la matriz A
        MPI_Scatter(
            A.data(), MATRIX_DIMENSION*rows_per_rank, MPI_DOUBLE, //datos de envio
            MPI_IN_PLACE, 0, MPI_DOUBLE, //datos de recepcion
            0,MPI_COMM_WORLD); // coordinador
        
        //--------enviar matriz B
        MPI_Bcast(B.data(),MATRIX_DIMENSION,MPI_DOUBLE,0,MPI_COMM_WORLD);

        //---------realizar el calculo C=A x b
        
        matrix_multi(A.data(),B.data(), C.data(), rows_per_rank, MATRIX_DIMENSION);

        //------recibir los resultados parciales

        MPI_Gather(MPI_IN_PLACE, 0, MPI_DOUBLE,
                    C.data(), rows_per_rank, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        C.resize(MATRIX_DIMENSION);

        std::printf("resultado\n");
        for(int i =0; i<MATRIX_DIMENSION;i++){
            std::printf("%.0f, ", C[i]);
        }


    }else{

        std::vector<double> A_local(MATRIX_DIMENSION*rows_per_rank);
        std::vector<double> B_local(MATRIX_DIMENSION);
        std::vector<double> C_local(rows_per_rank);

        //recibir matriz A
        MPI_Scatter(
            nullptr,0,MPI_DOUBLE,
            A_local.data(), MATRIX_DIMENSION*rows_per_rank, MPI_DOUBLE, //datos de recepcion
            0,MPI_COMM_WORLD); 

        std::printf("rank %d: [%.0f .. %.0f] \n",rank,A_local[0],A_local.back());

        //--------recibir matriz B
        MPI_Bcast(B_local.data(),MATRIX_DIMENSION,MPI_DOUBLE,0,MPI_COMM_WORLD);
        
        //---------realizar el calculo C=A x b
        int rows_per_rank_tmp = rows_per_rank;
        
        if(rank==nprocs-1){
            rows_per_rank_tmp = MATRIX_DIMENSION - rank * rows_per_rank;
        }

         matrix_multi(A_local.data(),B_local.data(), C_local.data(), rows_per_rank, MATRIX_DIMENSION);


        //enviar el resultado parcial
        MPI_Gather(C_local.data(), rows_per_rank, MPI_DOUBLE,
                    nullptr,0,MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    
    }

    MPI_Finalize();

    return 0;
}
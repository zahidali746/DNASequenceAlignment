/*-------------------------------------------------------------

 Copyright (C) 2000 Peter Clote. 
 All Rights Reserved.
 Copyright (C) 2010 Jawad Manzoor, Kiarash Rezahanjani, Erisa Dervishi. 

 Permission to use, copy, modify, and distribute this
 software and its documentation for NON-COMMERCIAL purposes
 and without fee is hereby granted provided that this
 copyright notice appears in all copies.


 THE AUTHOR AND PUBLISHER MAKE NO REPRESENTATIONS OR
 WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
 IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS
 AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED
 BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING
 THIS SOFTWARE OR ITS DERIVATIVES.

 -------------------------------------------------------------*/

/*************************************************

 Program for local sequence alignment, using the Smith-Waterman
 algorithm and assuming a LINEAR gap penalty.
 A traceback is used to determine the alignment, and
 is determined by choosing that direction, for
 which S(i-1,j-1)+sigma(a_i,b_j), S(i-1,j)+Delta and 
 S(i,j-1)+Delta is maximum, i.e.  for which 

 _
 |
 | H(i-1,j-1) + sigma(a_i,b_j)  (diagonal)
 H(i,j) =  MAX of   | H(i-1,j)   + delta           (up)
 | H(i,j-1)   + delta           (left)
 | 0
 | _


 is a maximum.

 *************************************************/

/*begin AMPP
 AMPP: SOME BUGS HAS BEEN SOLVED
 SOME CODE HAS BEEN IMPROVED 
 Code added or re-coded starts with begin AMPP and ends with
 end AMPP
 end AMPP 
 */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>   	// character handling
#include <stdlib.h>     // def of RAND_MAX 
/* begin AMPP */
/* Just a note:                       */
/* N must be the size of the arrays   */
/* Here is assume that the two arrays */
/* have the same size                 */

#define MAX_SEQ 50
#define PIPE_MSG 1
#define END_MSG 2

#define CHECK_NULL(_check) {\
   if ((_check)==NULL) \
      fprintf(stderr, "Null Pointer allocating memory\n");\
   }

#define AA 20           // number of amino acids
#define MAX2(x,y)     ((x)<(y) ? (y) : (x))
#define MAX3(x,y,z)   (MAX2(x,y)<(z) ? (z) : MAX2(x,y))

// function prototypes
void error(char *); /** error handling */

/* begin AMPP*/
int root, myId, mpiSize, a_b_size[2];
int char2AAmem[256];
int AA2charmem[AA];
int diag, down, right, max, Max, xMax, yMax, DELTA;
int *sendcnts; // Array of size mpiSize holding number of columns send to each node
int *displs; // Starting index 
short *a, *b, *a_partition; // protein a and b , and partition of protein a that every node will receive
int *hptr;
int **h;
int sim[AA][AA]; // PAM similarity matrix
double Ttot1, Ttot2, Tstart = 0, Tend = 0;
MPI_Status status;

void initChar2AATranslation(void);

Init(Argc,Argv)
int Argc; char **Argv;

{

	root=0;
	DELTA = atoi(Argv[4]);
	MPI_Init(&Argc,&Argv);

	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	MPI_Comm_rank(MPI_COMM_WORLD,&myId);

	if (myId == root)
	Ttot1 = MPI_Wtime(); //Start counting the clock.
}

main(argc,argv)
int argc; char **argv;
{
	Init(argc,argv);

	if(myId==root)
	{
		startMaster(argc,argv);
	}
	else
	{
		startSlave(argc,argv);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	fprintf(stderr, "END OF MAIN (rank: %d)\n", myId);
}

/* end AMPP*/

void startMaster(argc, argv)
int argc;
char **argv;
{
	// variable declarations
	FILE * in1, *in2, *pam;

	char ch;
	int temp;
	int i, j, tempi, tempj, x, y;
	int topskip, bottomskip;
	char *aout, *bout;
	int Aend, Bend, Abegin, Bbegin;

	// Max is first found maximum in similarity matrix H
	// max is auxilliary to determine largest of
	// diag,down,right, xMax,yMax are h-coord of Max

	//	short **xTraceback, **yTraceback;
	//	short *xTracebackptr, *yTracebackptr;
	int N;
	int nc;

	/**** Error handling for input file ****/
	if ((argc != 5) && (argc != 6)) {
		fprintf(stderr, "%s protein1 protein2 PAM gapPenalty [N]\n", argv[0]);
		exit(1);
	} else if (argc == 5) /* Maximum size of the proteins, they should have equal size */
	{
		/***** Initialization of input file and pair array **/
		in1 = fopen(argv[1], "r");
		in2 = fopen(argv[2], "r");
		pam = fopen(argv[3], "r");
		DELTA = atoi(argv[4]);
		N = MAX_SEQ;
	} else {
		in1 = fopen(argv[1], "r");
		in2 = fopen(argv[2], "r");
		pam = fopen(argv[3], "r");
		DELTA = atoi(argv[4]);
		N = atoi(argv[5]);
	}

	CHECK_NULL((aout = (char *) malloc(sizeof(char)*2*N)));
	CHECK_NULL((bout = (char *) malloc(sizeof(char)*2*N)));
	CHECK_NULL((a = (short *) malloc(sizeof(short)*(N+1))));
	CHECK_NULL((b = (short *) malloc(sizeof(short)*(N+1))));
	CHECK_NULL((hptr = (int *) malloc(sizeof(int)*(N+1)*(N+1))));
	CHECK_NULL((h = (int **) malloc(sizeof(int*)*(N+1))));
	/* Mount h[N][N] */
	for (i = 0; i <= N; i++)
	h[i] = hptr + i * (N + 1);

	initChar2AATranslation();

	/** read PAM250 similarity matrix **/
	fscanf(pam, "%*s");
	for (i = 0; i < AA; i++)
	for (j = 0; j <= i; j++) {
		if (fscanf(pam, "%d ", &temp) == EOF) {
			fprintf(stderr, "PAM file empty\n");
			fclose(pam);
			exit(1);
		}
		sim[i][j] = temp;
	}
	fclose(pam);
	for (i = 0; i < AA; i++)
	for (j = i + 1; j < AA; j++)
	{
		sim[i][j] = sim[j][i]; // symmetrify

	}

	/** read first file in array "a" **/
	i = 0;
	do {
		nc = fscanf(in1, "%c", &ch);
		if (nc > 0 && char2AAmem[ch] >= 0) {
			a[++i] = char2AAmem[ch];

		}
	}while (nc > 0 && (i < N));
	a[0] = i;
	fclose(in1);

	/** read second file in array "b" **/
	i = 0;
	do {
		nc = fscanf(in2, "%c", &ch);
		if (nc > 0 && char2AAmem[ch] >= 0) {
			b[++i] = char2AAmem[ch];
			// fprintf(stderr, "B %d - ", b[i]);
		}
	}while (nc > 0 && (i < N));
	b[0] = i;
	fclose(in2);

	/** initialize traceback array **/
	Max = xMax = yMax = 0;

	/* size of the proteings a and b | a will be partitioned | b the whole array is broadcasted */

	int rows, columns;

	a_b_size[0] = a[0];
	a_b_size[1] = b[0];

	// beginTime();
	MPI_Bcast(a_b_size, 2, MPI_INT, root, MPI_COMM_WORLD);
	//MPI_Barrier(MPI_COMM_WORLD);
	// endTime("Broad Cast sizes:");

	//beginTime();
	//bradcast size of the sequence a and b
	//---------------B------------------
	MPI_Bcast(b + 1, a_b_size[1], MPI_SHORT, root, MPI_COMM_WORLD);
	//endTime("Broad Cast B:");
	computeSendArray();
	CHECK_NULL((a_partition = (short *) malloc(sizeof(short)*(sendcnts[myId]))));

	//beginTime();
	//---------------A-------------------
	MPI_Scatterv(a, sendcnts, displs, MPI_SHORT, a_partition,
			sendcnts[myId], MPI_SHORT, root, MPI_COMM_WORLD);
	//endTime("Broad Cast A:");
	//beginTime();
	//-------------Score Matx--------------------
	MPI_Bcast(sim, AA * AA, MPI_INT, root, MPI_COMM_WORLD);
	//endTime("Broad Cast Score Mat:");

	int myColumns = getNoColumns(myId, mpiSize, a[0] + 1);

	for (i = 0; i <= a[0]; i++)
	h[i][0] = 0;
	for (j = 0; j <= b[0]; j++)
	h[0][j] = 0;

	for (i = 1; i <= b[0]; i++) {
		for (j = 1; j < myColumns; j++) {
			diag = h[i - 1][j - 1] + sim[b[i]][a[j]];
			down = h[i - 1][j] + DELTA;
			right = h[i][j - 1] + DELTA;
			max = MAX3(diag,down,right);
			if (max <= 0) {
				h[i][j] = 0;
			} else if (max == diag) {
				h[i][j] = diag;
			} else if (max == down) {
				h[i][j] = down;
			} else {
				h[i][j] = right;
			}
			if (max > Max) {
				Max = max;
				xMax = i;
				yMax = j;
			}
		} // end inner for loop

		MPI_Send( &h[i][j-1], 1, MPI_INT, myId + 1, PIPE_MSG, MPI_COMM_WORLD);
	}
	//Slave send the matrix back to the Master and Master combines the result into the final matrix.
	gatherFromSlaves(h, b[0]+1 );
}

void error(char * s) {
	fprintf(stderr, "%s\n", s);
	exit(1);
}

void initChar2AATranslation(void) {
	int i;
	for (i = 0; i < 256; i++)
		char2AAmem[i] = -1;
	char2AAmem['c'] = char2AAmem['C'] = 0;
	AA2charmem[0] = 'c';
	char2AAmem['g'] = char2AAmem['G'] = 1;
	AA2charmem[1] = 'g';
	char2AAmem['p'] = char2AAmem['P'] = 2;
	AA2charmem[2] = 'p';
	char2AAmem['s'] = char2AAmem['S'] = 3;
	AA2charmem[3] = 's';
	char2AAmem['a'] = char2AAmem['A'] = 4;
	AA2charmem[4] = 'a';
	char2AAmem['t'] = char2AAmem['T'] = 5;
	AA2charmem[5] = 't';
	char2AAmem['d'] = char2AAmem['D'] = 6;
	AA2charmem[6] = 'd';
	char2AAmem['e'] = char2AAmem['E'] = 7;
	AA2charmem[7] = 'e';
	char2AAmem['n'] = char2AAmem['N'] = 8;
	AA2charmem[8] = 'n';
	char2AAmem['q'] = char2AAmem['Q'] = 9;
	AA2charmem[9] = 'q';
	char2AAmem['h'] = char2AAmem['H'] = 10;
	AA2charmem[10] = 'h';
	char2AAmem['k'] = char2AAmem['K'] = 11;
	AA2charmem[11] = 'k';
	char2AAmem['r'] = char2AAmem['R'] = 12;
	AA2charmem[12] = 'r';
	char2AAmem['v'] = char2AAmem['V'] = 13;
	AA2charmem[13] = 'v';
	char2AAmem['m'] = char2AAmem['M'] = 14;
	AA2charmem[14] = 'm';
	char2AAmem['i'] = char2AAmem['I'] = 15;
	AA2charmem[15] = 'i';
	char2AAmem['l'] = char2AAmem['L'] = 16;
	AA2charmem[16] = 'l';
	char2AAmem['f'] = char2AAmem['F'] = 17;
	AA2charmem[17] = 'L';
	char2AAmem['y'] = char2AAmem['Y'] = 18;
	AA2charmem[18] = 'y';
	char2AAmem['w'] = char2AAmem['W'] = 19;
	AA2charmem[19] = 'w';
}

int getNoColumns(int myId, int mpiSize, int totalColumn) {
	int col = (totalColumn / mpiSize) + (totalColumn % mpiSize > myId);
	return col;
}

/* 
 Slave 

 */

void startSlave(argc, argv)
int argc;
char **argv;
{
	Ttot1= MPI_Wtime();
	double Tstart, Tend, Trec=0, Tsend=0, Tcomp=0, Tdata;
	//int a_b_size[2];//size of the proteings a and b | a will be partitioned | b the whole array is broadcasted
	int rows, columns=1, currentRecValue = 0, prevRecValue = 0;
	//MPI_Status status;

	MPI_Bcast(a_b_size, 2, MPI_INT, root, MPI_COMM_WORLD);

	//MPI_Barrier(MPI_COMM_WORLD);
	columns=getNoColumns(myId, mpiSize, a_b_size[0]+1);

	CHECK_NULL((a = (short *) malloc(sizeof(short))));
	CHECK_NULL((a_partition = (short *) malloc(sizeof(short) * (columns))));
	CHECK_NULL((b = (short *) malloc(sizeof(short) * (a_b_size[1]))));

	computeSendArray();

	//---------------B------------------
	MPI_Bcast(b, a_b_size[1], MPI_SHORT, root, MPI_COMM_WORLD);

	//---------------A-------------------
	MPI_Scatterv(a, sendcnts, displs, MPI_SHORT, a_partition,
			sendcnts[myId], MPI_SHORT, root, MPI_COMM_WORLD);

	//-------------Score Matx--------------------
	MPI_Bcast(sim, AA*AA, MPI_INT, root, MPI_COMM_WORLD);

	Tdata= MPI_Wtime();//time to receive data from master

	rows = a_b_size[1] + 1; //including 0 row in the beginning

	CHECK_NULL((hptr = (int *) malloc(sizeof(int) * (rows * columns))));
	CHECK_NULL((h = (int **) malloc(sizeof(int*) * (rows))));

	// Mount h[N][N] 
	int i, j;
	for (i = 0; i < rows; i++)
	h[i] = hptr + i * (columns);

	for (j = 0; j < columns; j++)
	h[0][j] = 0;
	//fprintf(stderr, "SLAVE Delta %d \n", DELTA);
	//fprintf(stderr, "\n checkpoint 1 Rank %d Rows: %d columns %d  \n",myId, rows,columns);
	//  Tstart= MPI_Wtime();
	for (i = 1; i < rows; i++) {
		Tstart= MPI_Wtime(); //receive time
		MPI_Recv(&currentRecValue, 1, MPI_INT, (myId-1), MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
		Tend= MPI_Wtime(); //end receive time
		Trec+=Tend-Tstart;

		Tstart= MPI_Wtime(); //computation time
		diag = prevRecValue + sim[b[i-1]][a_partition[0]];

		down = h[i - 1][0] + DELTA;
		right = currentRecValue + DELTA;
		max = MAX3(diag, down, right);
		//fprintf(stderr, "IND 0 down %d right %d diag %d max %d\n", down,right, diag, max);
		if (max <= 0) {
			h[i][0] = 0;

		} else if (max == diag) {
			h[i][0] = diag;

		} else if (max == down) {
			h[i][0] = down;

		} else {
			h[i][0] = right;

		}
		if (max > Max) {
			Max = max;
			xMax = i;
			yMax = 0;
		}
		prevRecValue = currentRecValue;

		for (j = 1; j < columns; j++) {

			diag = h[i - 1][j - 1] + sim[b[i-1]][a_partition[j]];
			down = h[i - 1][j] + DELTA;
			right = h[i][j - 1] + DELTA;

			max = MAX3(diag, down, right);
			//fprintf(stderr, "IND >0 down %d right %d diag %d max %d\n", down,right, diag, max);
			if (max <= 0) {
				h[i][j] = 0;

			} else if (max == diag) {
				h[i][j] = diag;

			} else if (max == down) {
				h[i][j] = down;

			} else {
				h[i][j] = right;

			}
			if (max > Max) {
				Max = max;
				xMax = i;
				yMax = j;
			}

		} // end inner for loop
		Tend= MPI_Wtime(); //computation tim 
		Tcomp+=Tend-Tstart;
		if (myId != mpiSize - 1)
		{
			//  fprintf(stderr, "\n Rank %d Sends Max %d H %d \n", myId, max,h[i][j-1]);                 
			Tstart= MPI_Wtime();
			MPI_Send(&max, 1, MPI_INT, myId + 1, PIPE_MSG, MPI_COMM_WORLD);
			Tend= MPI_Wtime();
			Tsend+=Tend-Tstart;

		}
	}//end outer for loop
	// Tend= MPI_Wtime();//Start counting the clock.


	Tstart= MPI_Wtime();
	MPI_Send(hptr, (rows) * (columns), MPI_INT, root, END_MSG, MPI_COMM_WORLD);
	Tend= MPI_Wtime();

	Ttot2= MPI_Wtime();

	printf("\n======= ==== Slave ==== =======\n");

	fprintf(stderr,"Slave: Rank: %d %s  %f  \n", myId, "Total Exec Time:", (float)(Ttot2-Ttot1));
	fprintf(stderr,"Slave: Rank: %d %s  %f  \n", myId, "Received data from Master:", (float)(Tdata-Ttot1));
	fprintf(stderr,"Slave: Rank: %d %s  %f %s %f \n", myId, "Avg Computation one row partition:", (float)(Tcomp/(rows-1)), "Total Computation", Tcomp);
	fprintf(stderr,"Slave: Rank: %d %s  %f  \n", myId, "Avg Time to receive from prev:", (float)(Trec/(rows-1)));
	fprintf(stderr,"Slave: Rank: %d %s  %f  \n", myId,"Avg Time to send to next:", (float)(Tsend/(rows-1)));
	fprintf(stderr,"Slave: Rank: %d %s  %f  \n", myId,"Time to send result to root:", (float)(Tend-Tstart));

}

void printArray(short *mat, int size, int rank) {
	int i;
	fprintf(stderr, "\n %d \n", rank);
	for (i = 0; i < size; i++) {
		fprintf(stderr, "%d ", mat[i]);

	}
	printf("\n \n");
}

void printMatrix(int **mat, int row, int col, int rank) {
	int i, j;
	fprintf(stderr, "\n %d \n", rank);

	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			fprintf(stderr, "%d ", mat[i][j]);

		}
		printf(" \n ");
	}
	printf("\n ");
}

void gatherFromSlaves(int **h, int rows) {

	double Tstart, Tend, Tps, Tpe, *waiting_time;
	//store wait time for each process to transfer data to master 
	CHECK_NULL((waiting_time= (double *) malloc(sizeof(double) * (mpiSize))));
	Tstart = MPI_Wtime();

	int rank, r, c, k, m, columns = sendcnts[0];//columns = columns assigned to root which certainly always has the max columns 
	rows = a_b_size[1] + 1;

	int *hptr_Part, **h_Part;
	CHECK_NULL((hptr_Part = (int *) malloc(sizeof(int) * (rows) * (columns))));
	CHECK_NULL((h_Part = (int **) malloc(sizeof(int*) * (rows))));

	for (rank = 1; rank < mpiSize; rank++) {
		Tps = MPI_Wtime();
		MPI_Recv(hptr_Part, (rows * sendcnts[rank]), MPI_INT, rank,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		Tpe = MPI_Wtime();
		waiting_time[rank] = Tpe - Tps;

		for (k = 0; k < rows; k++)
			h_Part[k] = hptr_Part + k * (sendcnts[rank]);

		for (r = 0; r < rows; r++)
			for (c = 0; c < sendcnts[rank]; c++) {

				int org_col = displs[rank] + c;
				h[r][org_col] = h_Part[r][c];

			}

	}
	Tend = MPI_Wtime();

	Ttot2 = MPI_Wtime();
	printf("\n======= Master: Gather Function =======\n");
	fprintf(stderr, "\nMaster:Total exe time of Master: %f  \n", (float) (Ttot2
			- Ttot1));
	fprintf(stderr, "Master:Total Gathering Time: %f  \n\n", (float) (Tend
			- Tstart));
	for (r = 1; r < mpiSize; r++) {
		fprintf(stderr,
				"Master:From Rank: %d -> Transfer time to master: %f \n", r,
				(float) waiting_time[r]);
	}

	printf("\n\n\n Final H Matrix \n");
	for (r = 0; r < rows; r++) {
		for (c = 0; c < rows; c++)
			fprintf(stderr, " %d ", h[r][c]);
		fprintf(stderr, "\n %s \n", "-");
	}
}

//compute base on size of a+1
void computeSendArray() {
	int i;

	CHECK_NULL((sendcnts = (int *) malloc(sizeof(int)*(mpiSize))));
	CHECK_NULL((displs = (int *) malloc(sizeof(int)*(mpiSize))));

	int offset = 0;

	for (i = 0; i < mpiSize; i++) {
		sendcnts[i] = getNoColumns(i, mpiSize, a_b_size[0] + 1);
		displs[i] = offset;
		offset += sendcnts[i];
	}
}

void beginTime() {
	Tstart = MPI_Wtime();
}
void endTime(char *msg) {
	Tend = MPI_Wtime();
	fprintf(stderr, "Master: %s  %f  \n", msg, (float) (Tend - Tstart));
}


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

//#define UNOPTIMIZED
#define OPTIMIZED

#define REG_A_SIZE 5
#define REG_B_SIZE 5
#define SNP_LENGTH 20

typedef unsigned int myDataType;
typedef double myResultType;

double StartTime;
double FinishTime;

double gettime(void)
{
	struct timeval ttime;
	gettimeofday(&ttime , NULL);
	return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

//initialize array with 2^16 elements
int bits_in_16bits [0x1u << 16];

//////////////////////////////////

void initialize_array(myDataType *array, int size, myDataType value)
{

    int count;
    // The count variable is used to parse through the array elements. In each for loop, count is incremented by 1.

    for (count=0; count<size; count++)
    {
        *(array+count) = value; //could also be written as array[count]
        //the name array is a pointer to the first element of the array
        //incrementing it gives us the pointer of the next element.
        //In truth the value of array doesnt increase by count but by 8*count. We are using the datatype long long unsigned so the distance from one array element to the next is 8 bytes of RAM.
    }

    //function of type void requires no return

}

//////////////////////////////////

/* LUT-based population counter */
int iterated_bitcount(unsigned int n)
/*iterated_bitcount:
	0x01 mask (0000 0001) is applied and the AND result is added to count
	n is then shifted by 1 to the right, which equals dividing n by 2^1
*/
{
	int count=0;

	while(n)
	{
		count += n & 0x1u ;
		n >>= 1 ;
	}

	return count;
}

void compute_bits_in_16bits(void)
/*compute_bits_in_16bits:
	fill in the empty array with the count of bits of each array address
	bits_in_16bits[1] = 1
	bits_in_16bits[2] = 2
	bits_in_16bits[3] = 2
	bits_in_16bits[4] = 3
	That will be our lookup table
*/
{
	unsigned int i;

	for (i = 0; i < (0x1u<<16); i++)
		bits_in_16bits[i] = iterated_bitcount(i);
}

int POPCOUNT (unsigned int n)
/*
bits_in_16bits [n         & 0xffffu] gives the bitcount of the first last 16 bits then after shifting by 16 we get the other 16 sum and add both together.
*/
{
	/* works only for 32-bit unsigned int*/
	return bits_in_16bits [n         & 0xffffu] //mask 0000 0000 0000 0000 1111 1111 1111 1111
	    +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}

#ifdef UNOPTIMIZED
myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b,  int size_a, int size_b, int size_snp, int snp_index_a, int snp_index_b)
{
	myResultType result = 0.0f;

	int i;
	int counter_1 = 0; // Counts 1s in snp a
	int counter_2 = 0; // Counts 1s in snp b
	int counter_3 = 0; // Counts pairs of 1s in the pair of snps

	for(i=0;i<size_snp;i++)
	{
		counter_1 += region_a[snp_index_a*size_snp+i]; //from given region select a specific snp and count its 1s
		counter_2 += region_b[snp_index_b*size_snp+i];
		counter_3 += region_a[snp_index_a*size_snp+i] & region_b[snp_index_b*size_snp+i];
	}

	if((counter_1==SNP_LENGTH)||(counter_2==SNP_LENGTH))
		return result;

	// Calculate the squared Pearson's correlation coefficient as measure of LD
	myResultType val_1 = ((myResultType)counter_1)/SNP_LENGTH;
	myResultType val_2 = ((myResultType)counter_2)/SNP_LENGTH;
	myResultType val_3 = ((myResultType)counter_3)/SNP_LENGTH;

	result = ((val_3-val_1*val_2)*(val_3-val_1*val_2));
	result /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));

	assert(result>=0.0000);
	assert(result<=1.0001);

	return result;
}

void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType * results, int size_a, int size_b, int size_snp)
{
	int i, j, k;//unnecessary k?
	for(i=0;i<size_a;i++)
	{
		for(j=0;j<size_b;j++)
		{
			//for all indices i and j compute LD
			results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_a, size_b, size_snp, i, j);
		}
	}
}
#endif

#ifdef OPTIMIZED
myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b,  int size_a, int size_b, int size_snp, int snp_index_a, int snp_index_b)
{
	myResultType result = 0.0f;

	// T2: LD computation on compressed data

	assert(result>=0.0000);
	assert(result<=1.0001);

	return result;
}

void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType * results, int size_a, int size_b, int size_snp)
{

	int i, j, k;

	// T1: In-place data compression
	//An expression a ? b : c evaluates to b if the value of a is true, and otherwise to c.

	int size_in_bits = sizeof(myDataType)*8; // 32 bits
	int sequence_size_in_compact_words = (size_snp/size_in_bits)+((size_snp%size_in_bits)!=0?1:0); //fit 20 into 32 bits-> only need 1. This is the size of the compressed array.

	//////////////////////////////////////////

	int totalsize_a = size_a*sequence_size_in_compact_words;
	myDataType compact_region_a[totalsize_a];
	initialize_array(compact_region_a, totalsize_a, 0); //initialize to 0

	int compact_word_index = 0; //variable used asindex to iterate the compressed arrays
	for(i=0;i<size_a;i++)
	{

		int bit_counter = 0;
		myDataType newBit = 0;
		myDataType tmpWord = 0;

		for(j=0;j<size_snp;j++)
		{
			newBit = (region_a[i*size_snp+j] & 1); //AND operation 1 is 00000...01
			tmpWord = tmpWord<<1; //Shift operation
			tmpWord = tmpWord | newBit;
			bit_counter++;

			if(bit_counter==size_in_bits)
			{
				compact_region_a[compact_word_index] = tmpWord;
				tmpWord = 0;
				bit_counter = 0;
				compact_word_index++;
			}

		}
		if(bit_counter!=size_in_bits)
			compact_region_a[compact_word_index] = tmpWord;
			compact_word_index++;
	}

	//////////////////////////////////////////

	myDataType * tmpmem_b = NULL;
	int totalsize_b = size_b*sequence_size_in_compact_words;
	myDataType compact_region_b[totalsize_b];

	compact_word_index = 0; //variable used asindex to iterate the compressed arrays
	for(i=0;i<size_b;i++)
	{
		int bit_counter = 0;
		myDataType newBit = 0;
		myDataType tmpWord = 0;

		for(j=0;j<size_snp;j++)
		{
			newBit = (region_b[i*size_snp+j] & 1); //AND operation 1 is 00000...01
			tmpWord = tmpWord<<1; //Shift operation
			tmpWord = tmpWord | newBit;
			bit_counter++;

			if(bit_counter==size_in_bits)
			{
				compact_region_b[compact_word_index] = tmpWord;
				tmpWord = 0;
				bit_counter = 0;
				compact_word_index++;
			}

		}
		if(bit_counter!=size_in_bits)
			compact_region_b[compact_word_index] = tmpWord;
			compact_word_index++;
	}
	printf("\n%u\n",compact_region_b[3]);
	return;
	//////////////////////////////////////////

	for(i=0;i<size_a;i++)
	{
		for(j=0;j<size_b;j++)
		{
			results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_a, size_b, size_snp, i, j);
		}
	}
}

#endif




int main()
{
	// Random seed for rand()
	int seed = 12345;
	srand(seed);

	// LUT-based population counter initialization
	compute_bits_in_16bits();

	// Memory allocation for two regions of SNPs
	myDataType * region_a = (myDataType*)malloc(sizeof(myDataType)*REG_A_SIZE*SNP_LENGTH);
	assert(region_a!=NULL);

	myDataType * region_b = (myDataType*)malloc(sizeof(myDataType)*REG_B_SIZE*SNP_LENGTH);
	assert(region_b!=NULL);

	int i;

	// Initialization with 0 or 1
	for(i=0;i<REG_A_SIZE*SNP_LENGTH;i++)
		region_a[i] = (myDataType)((rand()%2)); //rand %2 means binary, 0 or 1

	for(i=0;i<REG_B_SIZE*SNP_LENGTH;i++)
		region_b[i] = (myDataType)((rand()%2));

	int j;

	// Print content of region A
	printf("\nRegion A:\n");
	for(i=0;i<REG_A_SIZE;i++)
	{
		printf("SNP %d:\t", i);

		for(j=0;j<SNP_LENGTH;j++)
			printf("%u ", region_a[i*SNP_LENGTH+j]);

		printf("\n");

	}

	// Print content of region B
	printf("\nRegion B:\n");
	for(i=0;i<REG_B_SIZE;i++)
	{
		printf("SNP %d:\t", i);

		for(j=0;j<SNP_LENGTH;j++)
			printf("%u ", region_b[i*SNP_LENGTH+j]);

		printf("\n");

	}

	// Allocate memory for the results, i.e., LD scores
	myResultType * results = (myResultType*)malloc(sizeof(myResultType)*REG_B_SIZE*REG_A_SIZE);
	assert(results!=NULL);

	// Start timer
	StartTime = gettime();

	// Call function to compute all pairwise LD scores
	compute_all_pairwise_ld (region_a, region_b, results, REG_A_SIZE, REG_B_SIZE, SNP_LENGTH);
	return 1;

	// Stop timer
	FinishTime = gettime();

	printf("\nTotal execution time %.5f seconds\n", (FinishTime-StartTime));

	// Print LD scores and checksum
	myResultType results_checksum = 0.0f;

	printf("\nLD values\n");
	for(i=0;i<REG_A_SIZE;i++)
	{
		for(j=0;j<REG_B_SIZE;j++)
		{
			printf("(A%d, B%d) = %f\t", i, j, results[i*REG_B_SIZE+j]);
			results_checksum += results[i*REG_B_SIZE+j];
		}
		printf("\n");
	}

	printf("\nTotal LD = %.10f\n", results_checksum);

	// Free memory
	free(region_a);
	free(region_b);
	free(results);

	return 0;
}

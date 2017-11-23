#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

//#define UNOPTIMIZED

/////////////////////////

//#define OPTIMIZED
//#define OPTIMIZATION_REGISTER
#define OPTIMIZATION_REGISTER_LOOP_UNROLLING

/////////////////////////

#define REG_A_SIZE 1000
#define REG_B_SIZE 1000
#define SNP_LENGTH 10000

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

int bits_in_16bits [0x1u << 16];

/* LUT-based population counter */
int iterated_bitcount(unsigned int n)
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
{
	unsigned int i;

	for (i = 0; i < (0x1u<<16); i++)
	{
		bits_in_16bits[i] = iterated_bitcount(i);
	}
}

int POPCOUNT (unsigned int n)
{
	/* works only for 32-bit unsigned int*/
	return bits_in_16bits [n         & 0xffffu]
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
        counter_1 += region_a[snp_index_a*size_snp+i];
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
    int i, j;
    for(i=0;i<size_a;i++)
    {
        for(j=0;j<size_b;j++)
        {
            results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_a, size_b, size_snp, i, j);
        }
    }
}
#endif

#ifdef OPTIMIZED
void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType * results, int size_a, int size_b, int size_snp)
{

    int i, j;

    //variables for get pairwise scores:
    int compressed_size = (size_snp/(sizeof(myDataType)*8)+(size_snp%(sizeof(myDataType)*8)!=0?1:0));
    int k, counter_1, counter_2, counter_3;

    // (A) In-place Data Compression

    for(i=0;i<size_a;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_a[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_a[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_a[i*size_snp+newind] = temp;
    }

    for(i=0;i<size_b;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_b[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_b[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_b[i*size_snp+newind] = temp;
    }
    // End of In-place Data Compression



    // (B) Optimization target
    for(i=0;i<size_a;i++)
    {
        for(j=0;j<size_b;j++)
        {
            //results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_a, size_b, size_snp, i, j, compressed_size);

            counter_1 = counter_2 = counter_3 = 0;
            // bit-counting on compressed data
            for(k=0;k<compressed_size;k++)
            {//i = snp_index_a , j = snp_index_b
                counter_1 += POPCOUNT (region_a[i*size_snp+k]);
                counter_2 += POPCOUNT (region_b[j*size_snp+k]);
                counter_3 += POPCOUNT (region_a[i*size_snp+k] & region_b[j*size_snp+k]);
            }
            // End of LD computation on compressed data


            if((counter_1==SNP_LENGTH)||(counter_2==SNP_LENGTH))
                results[i*size_b+j]= 0.0f;

            myResultType val_1 = ((myResultType)counter_1)/SNP_LENGTH;
            myResultType val_2 = ((myResultType)counter_2)/SNP_LENGTH;
            myResultType val_3 = ((myResultType)counter_3)/SNP_LENGTH;

            results[i*size_b+j]  = ((val_3-val_1*val_2)*(val_3-val_1*val_2));
            results[i*size_b+j]  /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));


            assert(results[i*size_b+j] >=0.0000);
            assert(results[i*size_b+j] <=1.0001);

        }
    }
    // End of (B) Optimization target
}

#endif

#ifdef OPTIMIZATION_REGISTER
/*
The benefits of the register storage class are greatest for variables that the
function uses frequently, such as the counter variable for a loop.
The register keyword can be used only with simple numeric variables, not arrays or
structures. Also, it can’t be used with either static or external storage classes. You can’t define a pointer to a register variable.
*/

myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b, int size_snp, int i, int j)
{
    register int counter_1, counter_2, counter_3;
    counter_1 = counter_2 = counter_3 = 0;
    register int compressed_size = (size_snp/(sizeof(myDataType)*8)+(size_snp%(sizeof(myDataType)*8)!=0?1:0));
    register int k;
    myResultType result= 0.0f;

    // bit-counting on compressed data
    for(k=0;k<compressed_size;k++)
    {//i = snp_index_a , j = snp_index_b
        counter_1 += POPCOUNT (region_a[i*size_snp+k]);
        counter_2 += POPCOUNT (region_b[j*size_snp+k]);
        counter_3 += POPCOUNT (region_a[i*size_snp+k] & region_b[j*size_snp+k]);
    }
    // End of LD computation on compressed data

    if((counter_1==SNP_LENGTH)||(counter_2==SNP_LENGTH))
        return result;
    myResultType val_1 = ((myResultType)counter_1)/SNP_LENGTH;
    myResultType val_2 = ((myResultType)counter_2)/SNP_LENGTH;
    myResultType val_3 = ((myResultType)counter_3)/SNP_LENGTH;

    result  = ((val_3-val_1*val_2)*(val_3-val_1*val_2));
    result  /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));

    assert(result >=0.0000);
    assert(result <=1.0001);

    return result;
}



void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType * results, int size_a, int size_b, int size_snp)
{

    //For variables that are used to iterate loops it should be fruitful to have them on register
    register int i, j;

    // (A) In-place Data Compression

    for(i=0;i<size_a;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_a[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_a[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_a[i*size_snp+newind] = temp;
    }

    for(i=0;i<size_b;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_b[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_b[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_b[i*size_snp+newind] = temp;
    }
    // End of In-place Data Compression



    // (B) Optimization target
    for(i=0;i<size_a;i++)
    {
        for(j=0;j<size_b;j++)
        {
            results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j);

        }
    }
    // End of (B) Optimization target
}
#endif

#ifdef OPTIMIZATION_REGISTER_LOOP_UNROLLING
myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b, int size_snp, int i, int j)
{
    register int counter_1, counter_2, counter_3;
    counter_1 = counter_2 = counter_3 = 0;
    register int compressed_size = (size_snp/(sizeof(myDataType)*8)+(size_snp%(sizeof(myDataType)*8)!=0?1:0));
    register int k;
    myResultType result= 0.0f;

    // bit-counting on compressed data
    for(k=0;k<compressed_size;k++)
    {//i = snp_index_a , j = snp_index_b
        counter_1 += POPCOUNT (region_a[i*size_snp+k]);
        counter_2 += POPCOUNT (region_b[j*size_snp+k]);
        counter_3 += POPCOUNT (region_a[i*size_snp+k] & region_b[j*size_snp+k]);
    }
    // End of LD computation on compressed data

    if((counter_1==SNP_LENGTH)||(counter_2==SNP_LENGTH))
        return result;
    myResultType val_1 = ((myResultType)counter_1)/SNP_LENGTH;
    myResultType val_2 = ((myResultType)counter_2)/SNP_LENGTH;
    myResultType val_3 = ((myResultType)counter_3)/SNP_LENGTH;

    result  = ((val_3-val_1*val_2)*(val_3-val_1*val_2));
    result  /= (val_1*val_2*(1.0-val_1)*(1.0-val_2));

    assert(result >=0.0000);
    assert(result <=1.0001);

    return result;
}

void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType * results, int size_a, int size_b, int size_snp)
{

    //All of these variables are included in loops, which means it should be fruitful to have them on register
    register int i, j;

    // (A) In-place Data Compression

    for(i=0;i<size_a;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_a[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_a[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_a[i*size_snp+newind] = temp;
    }

    for(i=0;i<size_b;i++)
    {
        myDataType temp = 0;
        myDataType newbit = 0;
        int bitcounter = 0;
        int newind = 0;
        for(j=0;j<size_snp;j++)
        {
            newbit = region_b[i*size_snp+j]&1;
            temp <<= 1;
            temp |= newbit;
            bitcounter++;

            if(bitcounter==sizeof(myDataType)*8)
            {
                region_b[i*size_snp+newind] = temp;
                temp = 0;
                bitcounter=0;
                newind++;
            }

        }
        region_b[i*size_snp+newind] = temp;
    }
    // End of In-place Data Compression



    // (B) Optimization target
    int unroll_factor = 5;
    int unroll_stop = size_b-unroll_factor;
    int counter_1, counter_2, counter_3;
    myResultType val_1, val_2, val_3;
    for(i=0;i<size_a;i++)
    {
        for(j=0;j<unroll_stop;j+=unroll_factor)
        {

            results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j);
            results[i*size_b+j+1] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j+1);
            results[i*size_b+j+2] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j+2);
            results[i*size_b+j+3] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j+3);
            results[i*size_b+j+4] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j+4);
        }

        //finish residuals
        for(;j<size_b;j++)
        {
            results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j);

        }


    }
    // End of (B) Optimization target
}
#endif




#ifdef OPTIMIZATION_2_LOOP_UNROLLING_OPENMP
void compute_all_pairwise_ld (myDataType * region_a, myDataType * region_b, myResultType *
results, int size_a, int size_b, int size_snp)
{
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
		region_a[i] = (myDataType)((rand()%2));

	for(i=0;i<REG_B_SIZE*SNP_LENGTH;i++)
		region_b[i] = (myDataType)((rand()%2));

	int j;

/*	// Print content of region A
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
*/
	// Allocate memory for the results, i.e., LD scores
	myResultType * results = (myResultType*)malloc(sizeof(myResultType)*REG_B_SIZE*REG_A_SIZE);
	assert(results!=NULL);

	// Start timer
	StartTime = gettime();

	// Call function to compute all pairwise LD scores
	compute_all_pairwise_ld (region_a, region_b, results, REG_A_SIZE, REG_B_SIZE, SNP_LENGTH);

	// Stop timer
	FinishTime = gettime();

	printf("\nTotal execution time %.5f seconds\n", (FinishTime-StartTime));

	// Print LD scores and checksum
	myResultType results_checksum = 0.0f;

	//printf("\nLD values\n");
	for(i=0;i<REG_A_SIZE;i++)
	{
		for(j=0;j<REG_B_SIZE;j++)
		{
			//printf("(A%d, B%d) = %f\t", i, j, results[i*REG_B_SIZE+j]);
			results_checksum += results[i*REG_B_SIZE+j];
		}
		//printf("\n");
	}

	printf("\nTotal LD = %.10f\n", results_checksum);

	// Free memory
	free(region_a);
	free(region_b);
	free(results);

	return 0;
}

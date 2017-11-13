#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

//#define UNOPTIMIZED
#define OPTIMIZED

#define REG_A_SIZE 5
#define REG_B_SIZE 5
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
bits_in_16bits [n         & 0xffffu] gives the bitcount of the last 16 bits then after shifting by 16 we get the other 16 sum and add both together.
*/
{
    /* works only for 32-bit unsigned int*/
    return bits_in_16bits [n         & 0xffffu] //mask 0000 0000 0000 0000 1111 1111 1111 1111
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}

#ifdef UNOPTIMIZED
myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b, int size_snp, int snp_index_a, int snp_index_b)
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
            results[i*size_b+j] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j);
        }
    }
}
#endif

#ifdef OPTIMIZED
myResultType get_pairwise_ld_score (myDataType * region_a, myDataType * region_b, int size_snp, int snp_index_a, int snp_index_b, int compact_words_size)
{
    myResultType result = 0.0f;

    int i;
    int counter_1 = 0; // Counts 1s in snp a
    int counter_2 = 0; // Counts 1s in snp b
    int counter_3 = 0; // Counts pairs of 1s in the pair of snps

    // T2: LD computation on compressed data
    for(i=0;i<compact_words_size;i++)
    {
        counter_1 += POPCOUNT(region_a[snp_index_a+i]);
        counter_2 += POPCOUNT(region_b[snp_index_b+i]);
        counter_3 += POPCOUNT(region_a[snp_index_a+i]&region_b[snp_index_b+i]);

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

    //The benefits of the register storage class are greatest for variables that the function uses frequently, such as the counter variable for a loop.
    register int i, j;

    // T1: In-place data compression
    int size_in_bits = sizeof(myDataType)*8; // 32 bits
    int sequence_size_in_compact_words = (size_snp/size_in_bits)+((size_snp%size_in_bits)!=0?1:0); //Size of each compressed sequence. To fit 20 into 32 bits-> only need 1. To fit 100 into 32 -> 4

    //////////////////////////////////////////

    int totalsize_a = size_a*sequence_size_in_compact_words; //total size of compressed array
    int totalsize_b = size_b*sequence_size_in_compact_words;

    int compact_word_index = 0; //variable used asindex to iterate the compressed arrays

    // Loop Jamming: iterate both together until reaching the limit of the small one then iterate the bigger one alone.

    //But in order to do that the program needs to know which one is bigger, or if they are equal
    int big_size;
    int small_size;
    if (size_a > size_b)
    {
        big_size = size_a;
        small_size = size_b;
    }
    else if (size_a < size_b)
    {
        big_size = size_b;
        small_size = size_a;
    }
    else
    {
        small_size = big_size = size_a;
    }

    for(i=0;i<small_size;i++)
    {
        //
        int bit_counter = 0; // to count which bit we're on
        myDataType newBit_1 = 0; // to use when applying the mask
        myDataType newBit_2 = 0;
        myDataType tmpWord_1 = 0; // TmpWord is an intermediate variable for the compact form
        myDataType tmpWord_2 = 0;

        for(j=0;j<size_snp;j++)
        {
            newBit_1 = (region_a[i*size_snp+j] & 1); //AND operation 1 is 00000...01
            tmpWord_1 = tmpWord_1<<1; //Shift operation
            tmpWord_1 = tmpWord_1 | newBit_1; //OR operation

            newBit_2 = (region_b[i*size_snp+j] & 1); //AND operation 1 is 00000...01
            tmpWord_2 = tmpWord_2<<1; //Shift operation
            tmpWord_2 = tmpWord_2 | newBit_2; //OR operation

            bit_counter++;

            if(bit_counter==size_in_bits) //If bitcounter reaches this size we must use another spot in the array so we increment the index
            {
                region_a[compact_word_index] = tmpWord_1; //save the compact form
                region_b[compact_word_index] = tmpWord_2;
                tmpWord_1 = 0;
                tmpWord_2 = 0;
                bit_counter = 0;
                compact_word_index++;
            }


        }

        //If the bitcounter isn't equal to this size when it exits the loop (which is most likely) we need to add whatever is left.
        if(bit_counter!=size_in_bits)
            region_a[compact_word_index] = tmpWord_1;
            region_b[compact_word_index] = tmpWord_2;
            compact_word_index++;
    }

    // Repeat the process iterating for the remainder of snps of the bigger region
    if(size_a!=size_b)

        if(size_a>small_size)

            for(i=small_size;i<size_a;i++)
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
                        region_a[compact_word_index] = tmpWord;
                        tmpWord = 0;
                        bit_counter = 0;
                        compact_word_index++;
                    }

                }
                if(bit_counter!=size_in_bits)
                    region_a[compact_word_index] = tmpWord;
                    compact_word_index++;
            }

        else

            for(i=small_size;i<size_b;i++)
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
                        region_b[compact_word_index] = tmpWord;
                        tmpWord = 0;
                        bit_counter = 0;
                        compact_word_index++;
                    }

                }
                if(bit_counter!=size_in_bits)
                    region_b[compact_word_index] = tmpWord;
                    compact_word_index++;
            }
    //////////////////////////////////////////



    //////////////////////////////////////////

    register int counter_1;
    register int counter_2; //because we increment by sequence_size_in_compact_words, we also need 2 counters to be incremented by 1 to correctly address each result

    //If the SNP LENGTH is higher than 32 then it won't fit in 32 bits, in which case we need more than one spot in the array for each SNP. The number of spots needed is equal to the sequence_size_in_compact_words so i and j are incremented by the sequence_size_in_compact_words.
    counter_1=0;
    for(i=0;i<totalsize_a;i+=sequence_size_in_compact_words)
    {
        counter_2=0;
        for(j=0;j<totalsize_b;j+=sequence_size_in_compact_words)
        {
            results[counter_1*size_b+counter_2] = get_pairwise_ld_score (region_a, region_b, size_snp, i, j, sequence_size_in_compact_words);
            counter_2++;
        }
        counter_1++;
    }
    return;
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
    free(results);
    free(region_a);
    free(region_b);
    return 0;

}

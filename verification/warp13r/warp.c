/*********************************************
 * Reference implementation by WARP Team     *
**********************************************/
#include <stdio.h>
#include <stdlib.h>  // rand(), srand()
#include <time.h>    // time()
#include <sys/random.h>

// #define R 41    /*round number*/
#define RN    6    /*rotation number*/
#define BR    32  /*brunch number*/
#define BR_HALF    (BR / 2)  /*half of the branch number*/

int Sbox[BR_HALF] = { 0xc, 0xa, 0xd, 0x3, 0xe, 0xb, 0xf, 0x7, 0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6 };
int perm[BR] = { 31, 6, 29, 14, 1, 12, 21, 8, 27, 2, 3, 0, 25, 4, 23, 10, 15, 22, 13, 30, 17, 28, 5, 24, 11, 18, 19, 16, 9, 20, 7, 26, };
int RC0[41] = { 0x0U, 0x0U, 0x1U, 0x3U, 0x7U, 0xfU, 0xfU, 0xfU, 0xeU, 0xdU, 0xaU, 0x5U, 0xaU, 0x5U, 0xbU, 0x6U, 0xcU, 0x9U, 0x3U, 0x6U, 0xdU, 0xbU, 0x7U, 0xeU, 0xdU, 0xbU, 0x6U, 0xdU, 0xaU, 0x4U, 0x9U, 0x2U, 0x4U, 0x9U, 0x3U, 0x7U, 0xeU, 0xcU, 0x8U, 0x1U, 0x2U};
int RC1[41] = { 0x4U, 0xcU, 0xcU, 0xcU, 0xcU, 0xcU, 0x8U, 0x4U, 0x8U, 0x4U, 0x8U, 0x4U, 0xcU, 0x8U, 0x0U, 0x4U, 0xcU, 0x8U, 0x4U, 0xcU, 0xcU, 0x8U, 0x4U, 0xcU, 0x8U, 0x4U, 0x8U, 0x0U, 0x4U, 0x8U, 0x0U, 0x4U, 0xcU, 0xcU, 0x8U, 0x0U, 0x0U, 0x4U, 0x8U, 0x4U, 0xcU};

void enc(int *m, int *c, int *k, int nr);
void sboxkey(int *state, int *k, int r);
void permutation(int *state);

#define PRINT_INTER 0

void init_prng() {
    unsigned int initial_seed = -1;
    ssize_t temp;
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
    srand(initial_seed);   // Initialization, should only be called once. int r = rand();
	printf("[x] PRNG initialized by %lu random bytes: 0x%08X\n", temp, initial_seed);
}

void print_target_bits(int *collected_ciphers, int ns, int *selected_positions)
{

    int random_position;
    int nibble;
    int bp;
    int target_bit;
    printf("target positions: ");
    for (int bit = 0; bit < ns; bit++)
    {
        nibble = selected_positions[bit] / 4;
        bp = selected_positions[bit] % 4;
        target_bit = ((collected_ciphers[nibble]) >> (3 - bp)) & 0x1;
        printf("%d", target_bit);
    }
    printf("\nrandom positions: ");
    for (int bit = 0; bit < ns; bit++)
    {
        random_position = rand() % 128;
        nibble = random_position / 4;
        bp = random_position % 4;
        target_bit = ((collected_ciphers[nibble]) >> (4 - bp)) & 0x1;
        printf("%d", target_bit);
    }
    printf("\n");
}

void print_state(int *m)
{
    for (int i = 0; i < BR; i++)
    {
        printf("%x ", m[i]&0xf);
    }
    printf("\n");
};

void printState(int *state)
{
    printf("L: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 0]);
    }
    printf("R: ");
    for (int x = 0; x < BR_HALF; x++)
    {
        printf("%x ", state[2 * x + 1]);
    }
    printf("\n");
}

int main()
{
    init_prng();
    // Randomly generate master key
    int k[BR];
    for (int i = 0; i < BR; i++)
        k[i] = rand() & 0xf;
    // Randomly generate a plaintext
    int m[BR];
    for (int i = 0; i < BR; i++)
        m[i] = rand() & 0xf;
    // Randomly generate a ciphertext
    int c[BR];
    for (int i = 0; i < BR; i++)
        c[i] = rand() & 0xf;
    int collected_ciphers[BR] = {0};

    //####################################################
    //####################################################
    //####################################################
    // Output of our Python program:
    int NUMBER_OF_ROUNDS = 13;
    int na = 3;
    int nb = 20;
    int active_nibbles[3] = {4, 5, 11};
    int balanced_positions[20] = {4, 5, 6, 7, 36, 37, 38, 39, 60, 61, 62, 63, 116, 117, 118, 119, 124, 125, 126, 127};
    //####################################################
    //####################################################
    //####################################################

    // /*key*/
    // printf("key :\n");
    // print_state(k);
    // /*plaintext*/
    // printf("plaintext :\n");
    // print_state(m);
    // enc(m, c, k, NUMBER_OF_ROUNDS);
    // /*ciphertext*/
    // printf("ciphertext :\n");
    // print_state(c);
    // printf("\n");
    unsigned long long int upper_bound = 1 << (4*na);
    for (unsigned long long int ti = 0; ti < upper_bound; ti++)
    {
        for (int nibble = 0; nibble < na; nibble++)
        {
            m[active_nibbles[nibble]] = (ti >> 4*nibble) & 0xf;
        }
        enc(m, c, k, NUMBER_OF_ROUNDS);
        // print_state(c);
        for (int output_nibble = 0; output_nibble < BR; output_nibble++)
        {               
            collected_ciphers[output_nibble] ^= c[output_nibble];            
        }
    }
    print_state(collected_ciphers);
    print_target_bits(collected_ciphers, nb, balanced_positions);
    return 0;
}

void enc(int *m, int *c, int *k, int R)
{
    /*intermediate value*/
    int state[BR];

    /*left half intermediate value*/
    int temp[BR_HALF];

    for (int i = 0; i < BR; i++)
    {
        state[i] = m[i];
    }

    /*round function(1 to 40 round)*/
    for (int i = 0; i < R; i++)
    {
        #if PRINT_INTER
        printf("%d round\n", i + 1);
        printState(state);
        #endif

        for (int j = 0; j < BR_HALF; j++)
        {
            temp[j] = state[j * 2];
        }
        /*insert key and Sbox*/
        sboxkey(temp, k, i);
        /*XOR*/
        for (int j = 0; j < BR_HALF; j++)
        {
            state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
        }
        /*add round constants*/
        state[1] = state[1] ^ RC0[i];
        state[3] = state[3] ^ RC1[i];

        /*permutation*/
        permutation(state);
    }

    /*last round function */
    #if PRINT_INTER
    printf("%d round\n", R);
    printState(state);
    #endif
    
    // Deactivate last partial round
    // for (int j = 0; j < BR_HALF; j++)
    // {
    //     temp[j] = state[j * 2];
    // }

    // /*input key and  Sbox*/
    // sboxkey(temp, k, R - 1);
    // /*xor*/
    // for (int j = 0; j < BR_HALF; j++)
    // {
    //     state[2 * j + 1] = state[2 * j + 1] ^ temp[j];
    // }
    // /*add round constants*/
    // state[1] = state[1] ^ RC0[R-1];
    // state[3] = state[3] ^ RC1[R-1];

    #if PRINT_INTER
    printState(state);
    #endif

    /*no permutation in the last round*/

    /*copy ciphertext*/
    for (int i = 0; i < BR; i++)
    {
        c[i] = state[i];
    }

}

void sboxkey(int *state, int *k, int r)
{
    for (int i = 0; i < BR_HALF; i++)
    {
        state[i] = Sbox[state[i]] ^ k[(r % 2) * 16 + i];
    }
}

void permutation(int *state)
{
    int temp[BR];
    for (int j = 0; j < BR; j++)
    {
        temp[j] = state[j];
    }
    for (int j = 0; j < BR; j++)
    {
        state[perm[j]] = temp[j];
    }
}

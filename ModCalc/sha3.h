//https://github.com/jedisct1/TweetFIPS202
#ifndef SHA3_H

#ifdef __cplusplus
extern "C"
{
#endif
	
void FIPS202_SHA3_224(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);
void FIPS202_SHA3_256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);
void FIPS202_SHA3_384(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);
void FIPS202_SHA3_512(const unsigned char *input, unsigned int inputByteLen, unsigned char *output);
void FIPS202_SHAKE128(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);
void FIPS202_SHAKE256(const unsigned char *input, unsigned int inputByteLen, unsigned char *output, int outputByteLen);

#ifdef __cplusplus
}
#endif

#endif
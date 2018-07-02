typedef void MT19937;

#ifdef __cplusplus
extern "C" {
#endif

	MT19937 *Mt19937Alloc(int i);
	unsigned int Mt19937NextInt(const MT19937 *t);
	double Mt19937NextDouble(const MT19937 *t);
	void Mt19937Free(MT19937 *t);
#ifdef __cplusplus
}
#endif

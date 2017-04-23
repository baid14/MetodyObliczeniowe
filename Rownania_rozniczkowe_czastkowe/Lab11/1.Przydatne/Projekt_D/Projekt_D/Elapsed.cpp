#include "Elapsed.h"

CElapsed::CElapsed()
{
	// get the frequency of the counter
	m_iInitialized = QueryPerformanceFrequency( (LARGE_INTEGER *)&m_iFrequency );
}

CElapsed::~CElapsed()
{

}

int CElapsed::Begin()    // start timing
{
	if (!m_iInitialized)
		return 0;   // error - couldn't get frequency

	// get the starting counter value
	return QueryPerformanceCounter( (LARGE_INTEGER *)&m_iBeginTime );
}

double CElapsed::End()    // stop timing and get elapsed time in miliseconds
{

	// get the ending counter value
	QueryPerformanceCounter( (LARGE_INTEGER *)&endtime );

	// determine the elapsed counts
	__int64 elapsed = endtime - m_iBeginTime;

	// convert counts to time in seconds and return it
	m_dElapsed = (double)elapsed / (double)m_iFrequency;
	return m_dElapsed;
}

int CElapsed::Available()  // returns true if the perf counter is available
{ 
	return m_iInitialized; 
}

__int64 CElapsed::GetFreq() // return perf counter frequency as large int
{
	return m_iFrequency; 
}
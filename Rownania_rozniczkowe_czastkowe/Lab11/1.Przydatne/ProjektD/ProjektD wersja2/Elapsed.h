#pragma once
#include "stdafx.h"

class CElapsed
{
public:
	CElapsed();
	virtual ~CElapsed();

	int Begin();    // start timing
	double End();
	int Available();
	__int64 GetFreq();

	double m_dElapsed;
private :
    int m_iInitialized;
    __int64 m_iFrequency;
    __int64 m_iBeginTime;
	__int64 endtime;

};
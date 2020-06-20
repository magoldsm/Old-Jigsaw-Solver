#pragma once

class CProgressBar
{
public:
	CProgressBar() : m_Percent(0.0), m_Bars(nullptr) {}
	CProgressBar(int n) { m_Bars = new CProgressBar[n]; }
	~CProgressBar() { if (m_Bars) delete m_Bars; }

	void MakeSubbars(int n) { m_Bars = new CProgressBar[n]; }

	double			m_Percent;
	std::string		m_Text;
	LARGE_INTEGER	m_Time;
	CProgressBar*	m_Bars;
};


class CProgress
{
public:
	CProgress() : m_Bars(nullptr) {}
	CProgress(int n) { m_Bars = new CProgressBar[n]; }
	~CProgress() { if (m_Bars) delete m_Bars; }

	CProgressBar*	m_Bars;

	CProgressBar& operator[](int n) { return m_Bars[n]; }
	CProgressBar& operator()(int n) { return m_Bars[n]; }
	CProgressBar& operator()(int n, int m) { return m_Bars[n].m_Bars[m]; }

	void RestartReport(int n, bool bStart);				// Set bStart=1 to start a bar, 0 to reset it
	void RestartReport(int n, int m, bool bStart);		// Set bStart=1 to start a bar, 0 to reset it
	void UpdateReport();								// Update bars 
};

#define	PROGRESS_EUCLID		0
#define	PROGRESS_BIVERTEX	1
#define	PROGRESS_PLACING	2
#define	PROGRESS_COMPARING	3,0
#define	PROGRESS_CHECKING	3,1

extern CProgress Progress;
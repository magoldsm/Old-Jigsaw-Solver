
// Jigsaw Solver WDlg.h : header file
//

#pragma once
#include "CParameters.h"
#include "CPlot.h"

class CProgressBar;

// CJigsawSolverWDlg dialog
class CJigsawSolverWDlg : public CDialogEx, public CParameters
{
// Construction
public:
	CJigsawSolverWDlg(CWnd* pParent = nullptr);	// standard constructor

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_JIGSAWSOLVERW_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

public:
	CString m_PuzzleFile;
	CProgressCtrl m_barEuclidean;
	CProgressCtrl m_barBivertex;
	CProgressCtrl m_barPlacing;
	CProgressCtrl m_barComparing;
	CProgressCtrl m_barChecking;


	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:

public:
	afx_msg void OnEnKillfocusJstar();
	afx_msg void OnBnClickedSolve();
	afx_msg LRESULT OnProgress(WPARAM wParam, LPARAM lParam);
	afx_msg LRESULT OnErase(WPARAM wParam, LPARAM lParam);
	afx_msg LRESULT OnPlot(WPARAM wParam, LPARAM lParam);
	afx_msg LRESULT OnDelete(WPARAM wParam, LPARAM lParam);

	virtual BOOL OnCommand(WPARAM wParam, LPARAM lParam);
	void UpdateParameters();

	void UpdateBar(CProgressCtrl& ctrl, CProgressBar & bar);

	CPlot m_Plot;
	BOOL m_bShowPScores;
};

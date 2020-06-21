
// Jigsaw Solver WDlg.cpp : implementation file
//

#include "pch.h"
#include "Jigsaw Solver W.h"
#include "Jigsaw Solver WDlg.h"
#include "afxdialogex.h"
#include "CProgress.h"
#include "CParameters.h"
#include "Savitsy-Golay.h"
#include "Euclidean Signature.h"
#include "BivertexArcDecomposition.h"
#include "Utilities.h"
#include "SignatureSimilarity.h"

using namespace std;
using namespace cv;
using namespace Eigen;


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

std::vector<CPiece> Pieces;
Eigen::VectorXd Weights;
Eigen::VectorXd smoothVec, d1Vec, d2Vec;
double Dx, Dy, Dkappa, Dkappas;
double AverageLength;
long AverageSize;
CProgress Progress(N_PROGRESS);

void CProgress::UpdateReport()
{
	::PostMessage(m_hWndGUI, WM_PROGRESS, 0, 0L);
}


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CJigsawSolverWDlg dialog



CJigsawSolverWDlg::CJigsawSolverWDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_JIGSAWSOLVERW_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CJigsawSolverWDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SG_ORDER, m_nSGOrder);
	DDX_Text(pDX, IDC_SG_WINDOW, m_nSGWindow);
	DDX_Text(pDX, IDC_ALPHA, m_nAlpha);
	DDX_Text(pDX, IDC_BETA, m_nBeta);
	DDX_Text(pDX, IDC_GAMMA, m_nGamma);
	DDX_Text(pDX, IDC_C1, m_nC1);
	DDX_Text(pDX, IDC_C2, m_nC2);
	DDX_Text(pDX, IDC_K1, m_dK1);
	DDX_Text(pDX, IDC_K2, m_dK2);
	DDX_Text(pDX, IDC_K4, m_dK4);
	DDX_Text(pDX, IDC_LAMBDA0, m_nLambda0);
	DDX_Text(pDX, IDC_LAMBDA1, m_nLambda1);
	DDX_Text(pDX, IDC_NU, m_nNu);
	DDX_Text(pDX, IDC_EPSILON, m_dEpsilon);
	DDX_Text(pDX, IDC_RHO, m_dRho);
	DDX_Text(pDX, IDC_JMAX, m_nJmax);
	DDX_Text(pDX, IDC_ETA1, m_dEta1);
	DDX_Text(pDX, IDC_ETA2, m_dEta2);
	DDX_Text(pDX, IDC_Q1, m_dQ1);
	DDX_Text(pDX, IDC_Q2, m_dQ2);
	DDX_Text(pDX, IDC_Q2STAR, m_dQ2Star);
	DDX_Text(pDX, IDC_Q3, m_dQ3);
	DDX_Text(pDX, IDC_P0, m_P0[0]);
	DDX_Text(pDX, IDC_M0, m_M0[0]);
	DDX_Text(pDX, IDC_MU0, m_MU0[0]);
	DDX_Text(pDX, IDC_K30, m_K3[0]);
	DDX_Text(pDX, IDC_P1, m_P0[1]);
	DDX_Text(pDX, IDC_M1, m_M0[1]);
	DDX_Text(pDX, IDC_MU1, m_MU0[1]);
	DDX_Text(pDX, IDC_K31, m_K3[1]);
	DDX_Text(pDX, IDC_P2, m_P0[2]);
	DDX_Text(pDX, IDC_M2, m_M0[2]);
	DDX_Text(pDX, IDC_MU2, m_MU0[2]);
	DDX_Text(pDX, IDC_K32, m_K3[2]);
	DDX_Text(pDX, IDC_P3, m_P0[3]);
	DDX_Text(pDX, IDC_M3, m_M0[3]);
	DDX_Text(pDX, IDC_MU3, m_MU0[3]);
	DDX_Text(pDX, IDC_K33, m_K3[3]);
	DDX_Text(pDX, IDC_P4, m_P0[4]);
	DDX_Text(pDX, IDC_M4, m_M0[4]);
	DDX_Text(pDX, IDC_MU4, m_MU0[4]);
	DDX_Text(pDX, IDC_K34, m_K3[4]);
	DDX_Text(pDX, IDC_P5, m_P0[5]);
	DDX_Text(pDX, IDC_M5, m_M0[5]);
	DDX_Text(pDX, IDC_MU5, m_MU0[5]);
	DDX_Text(pDX, IDC_K35, m_K3[5]);
	DDX_Text(pDX, IDC_P6, m_P0[6]);
	DDX_Text(pDX, IDC_M6, m_M0[6]);
	DDX_Text(pDX, IDC_MU6, m_MU0[6]);
	DDX_Text(pDX, IDC_K36, m_K3[6]);
	DDX_Text(pDX, IDC_P7, m_P0[7]);
	DDX_Text(pDX, IDC_M7, m_M0[7]);
	DDX_Text(pDX, IDC_MU7, m_MU0[7]);
	DDX_Text(pDX, IDC_K37, m_K3[7]);
	DDX_Text(pDX, IDC_JSTAR, m_nJStar);
	if (!pDX->m_bSaveAndValidate)
		m_PuzzleFile = m_szPath;
	DDX_Text(pDX, IDC_PUZZLE_FILE, m_PuzzleFile);
	if (pDX->m_bSaveAndValidate)
		strcpy_s(m_szPath, sizeof(m_szPath), CT2A(m_PuzzleFile));
	DDX_Control(pDX, IDC_EUCLIDEAN_SIGNATURES, m_barEuclidean);
	DDX_Control(pDX, IDC_BIVERTEX_ARC, m_barBivertex);
	DDX_Control(pDX, IDC_PLACING_PIECES, m_barPlacing);
	DDX_Control(pDX, IDC_COMPARING_PIECES, m_barComparing);
	DDX_Control(pDX, IDC_CHECKING_FITS, m_barChecking);
}

BEGIN_MESSAGE_MAP(CJigsawSolverWDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_EN_KILLFOCUS(IDC_JSTAR, &CJigsawSolverWDlg::OnEnKillfocusJstar)
	ON_BN_CLICKED(IDC_SOLVE, &CJigsawSolverWDlg::OnBnClickedSolve)
	ON_MESSAGE(WM_PROGRESS, &CJigsawSolverWDlg::OnProgress)
END_MESSAGE_MAP()


// CJigsawSolverWDlg message handlers

BOOL CJigsawSolverWDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	Progress(3).MakeSubbars(2);
	Progress.m_hWndGUI = m_hWnd;	// Have the progress bar send us messages

	m_barEuclidean.SetRange(0, 1000);
	m_barBivertex.SetRange(0, 1000);
	m_barPlacing.SetRange(0, 1000);
	m_barComparing.SetRange(0, 1000);
	m_barChecking.SetRange(0, 1000);

	OnEnKillfocusJstar();

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CJigsawSolverWDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CJigsawSolverWDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CJigsawSolverWDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

#define	CTRLS(n)	IDC_P##n, IDC_M##n, IDC_MU##n, IDC_K3##n

int Ctrls[][4] = {
	{ CTRLS(0) },
	{ CTRLS(1) },
	{ CTRLS(2) },
	{ CTRLS(3) },
	{ CTRLS(4) },
	{ CTRLS(5) },
	{ CTRLS(6) },
	{ CTRLS(7) },
};

void CJigsawSolverWDlg::OnEnKillfocusJstar()
{
	UpdateData(TRUE);

	for (int i = 0; i < 8; i++)
	{
		if (i < m_nJStar)
		{
			GetDlgItem(Ctrls[i][0])->ShowWindow(SW_SHOW);
			GetDlgItem(Ctrls[i][1])->ShowWindow(SW_SHOW);
			GetDlgItem(Ctrls[i][2])->ShowWindow(SW_SHOW);
			GetDlgItem(Ctrls[i][3])->ShowWindow(SW_SHOW);
		}
		else
		{
			GetDlgItem(Ctrls[i][0])->ShowWindow(SW_HIDE);
			GetDlgItem(Ctrls[i][1])->ShowWindow(SW_HIDE);
			GetDlgItem(Ctrls[i][2])->ShowWindow(SW_HIDE);
			GetDlgItem(Ctrls[i][3])->ShowWindow(SW_HIDE);
		}
	}
}


BOOL CJigsawSolverWDlg::OnCommand(WPARAM wParam, LPARAM lParam)
{
	UINT notificationCode = (UINT)HIWORD(wParam);

	// For List control I handle it in another way....
	if ((notificationCode == EN_KILLFOCUS) ||
		(notificationCode == LBN_KILLFOCUS) ||
		(notificationCode == CBN_KILLFOCUS) ||
		(notificationCode == NM_KILLFOCUS) ||
		(notificationCode == WM_KILLFOCUS)) {

		CWnd *pFocus = CWnd::GetFocus();
		// call to a static function 

		// If we are changing the focus to another
		// control of the same window... 

		if (pFocus && (pFocus->GetParent() == this))
		{
			// Ok, if the focus is not in the cancel button...
			if (pFocus->GetDlgCtrlID() != IDCANCEL) {
				UpdateData(TRUE);
				UpdateParameters();
			}
		}
	}

	return CDialogEx::OnCommand(wParam, lParam);
}


void CJigsawSolverWDlg::UpdateParameters()
{
	SaveParams();
}


void CJigsawSolverWDlg::OnBnClickedSolve()
{
	::AfxBeginThread([](LPVOID pParam)->UINT {
		Solver();
		return 0;
	}, nullptr);
}

LRESULT CJigsawSolverWDlg::OnProgress(WPARAM wParam, LPARAM lParam)
{
	UpdateBar(m_barEuclidean, Progress(PROGRESS_EUCLID));
	UpdateBar(m_barBivertex, Progress(PROGRESS_BIVERTEX));
	UpdateBar(m_barPlacing, Progress(PROGRESS_PLACING));
	UpdateBar(m_barComparing, Progress(PROGRESS_COMPARING));
	UpdateBar(m_barChecking, Progress(PROGRESS_CHECKING));

	/*
	
	    if(isempty(timerhs))
        tstr = '';
    else
        if(ps == 1)
            if(isa(timerhs, 'uint64'))
                timerhs = toc(timerhs);
            end            
            tstr = ['Ran ' sec2str(timerhs)];           
        else
            if(ps == 0)
                ttime = toc(timerhs);           
                tstr = {[sec2str(ttime) ' Running' ]};
            else
                ttime = toc(timerhs);           
                tstr = {[sec2str(ttime) ' Running' ]; ['~' sec2str(ttime/ps-ttime) ' Remaining']};
            end
        end
    end
    set(ths(1, 1), 'String', [ num2str(floor(1000*ps)/10) '% Complete' ]);
    set(ths(1, 2), 'String', tstr);
    set(ahs, 'Children', [ths(1, 1) ths(1, 2) shs]);

	*/

	return 0;
}


void CJigsawSolverWDlg::UpdateBar(CProgressCtrl& ctrl, CProgressBar & bar)
{
	int nLower, nUpper;
	ctrl.GetRange(nLower, nUpper);
	int nDiff = nUpper - nLower;
	int val = (bar.m_Percent * nDiff) + nLower;

	ctrl.SetPos(val);
	if (bar.m_Percent != 0)
		TRACE("%d\n", val);
}


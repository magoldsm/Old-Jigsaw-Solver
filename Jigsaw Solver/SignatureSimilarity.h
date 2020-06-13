#pragma once

class CTracker
{
public:
	CTracker() {}

	void SetSize(size_t sz)
	{
		m_PlacedPieces.resize(sz);
		m_RemainingPieces.resize(sz);
		m_InactiveArcs.resize(sz);
		m_ActiveArcs.resize(sz);
		m_ActivePoints.resize(sz);
		m_Pc2Place.resize(sz);
	}
	std::vector<int>				m_PlacedPieces;
	std::vector<int>				m_RemainingPieces;
	std::vector<Eigen::VectorXi>	m_InactiveArcs;
	std::vector<Eigen::VectorXi>	m_ActiveArcs;
	std::vector<Eigen::VectorXi>	m_ActivePoints;
	std::vector<int>				m_Pc2Place;
	Curve							m_SolvedPuzzleBoundary;
	//std::vector<xxx>		m_InactivePoints;

};

struct GTransform
{
	GTransform() { theta = 0.0; dx = 0.0; dy = 0.0; }
	double		theta;
	double		dx;
	double		dy;
};


class CFit
{
public:
	CFit() : m_Size(0), m_Slot(0), m_Score(-1) {}

	Eigen::Vector2i					m_Pieces;
	int								m_Size;
	Eigen::MatrixX2i				m_Arcs;				// Rows given by m_Size
	std::vector<GTransform>			m_ArcTrans;			// Rows given by m_Size
	GTransform						m_gFit;
	double							m_Score;
	int								m_Slot;
};


class CPlacement
{
public:
	CPlacement(int nPiece, const Eigen::Vector4d& Score, const GTransform& gLock, const CFit& Fit)
		: m_nPiece(nPiece)
		, m_Score(std::move(Score))
		, m_gLock(std::move(gLock))
		, m_Fit(std::move(Fit))
	{}
	int								m_nPiece;
	Eigen::Vector4d					m_Score;
	GTransform							m_gLock;
	CFit							m_Fit;
	//***BUGBUG***					m_Neighbors;
};

void PlacePieces();


class CPScore
{
public:
	CPScore(size_t sz);
	CPScore() : m_Size(0), m_ArcScores(NULL) {}
	~CPScore();

	void SetSize(size_t sz);

	bool IsEmpty(size_t nRow, size_t nCol) { return m_ArcScores[nRow*m_Size + nCol].size() == 0; }

	Eigen::MatrixXd& operator()(size_t nRow, size_t nCol)
	{
		return (m_ArcScores[nRow*m_Size + nCol]);
	}

	void Display();
	void Display(size_t nRow, size_t nCol);

private:
	size_t				m_Size;
	Eigen::MatrixXd*	m_ArcScores;
};
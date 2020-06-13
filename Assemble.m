% References: 
% 1. Hoff, D., Olver, P.J., Automatic Solution of Jigsaw Puzzles, preprint, University of Minnesota, 2012.
% 2. Hoff, D., Olver, P.J.: Extensions of Invariant Signatures for Object Recognition. J. Math. Imaging Vis., Online First(TM), 7 June 2012.


% This function attemps to solve a puzzle using the algorithm 
% described in [1]. It is intended for use by the Solve_Puzzle() 
% function
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'puzzle':   This should be a single row cell array, each of whose 
            entries represents the boundary of a puzzle piece. Each 
            boundary should be a discretized planar curve of m points 
            represented as an m-by-2 matrix, each row of which 
            specifies a point in R^2. The curves are assumed to be 
            closed; no repetition of initial points is necessary. 

'plotter':  This input should belong to the set {0, 1, 2, 3} and 
            controls what aspects of the algorithm are visualized as 
            the code runs. In all cases, a progress report is 
            generated. If plotter = 0, no other figures are 
            generated. If plotter = 1, the final solution is also 
            plotted. If plotter = 2, the current subpuzzle is plotted 
            as the puzzle is assemble. If plotter = 3, each 
            piece-locking is also plotted.  

'saver':    This logical input determines whether or not data is
            saved as the algorithm runs. Data is saved to the current
            directory as 'CurrentPuzzleData.mat' and can be used with
            the 'Assemble.m' function to resume computation at the 
            point of saving.

'parameters':   This variable contains the parameters (described in 
            [1]) used by the algorithm. Use the Solve_Puzzle() 
            function for a GUI to input these parameters, or see its 
            code for the form this input must take.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'pieces':   This output is a 1-by-n structure each of whose columns
            corresponds to a puzzle piece. The 'Points' field 
            contains the boundary of a piece as a planar curve of m 
            points represented as an m-by-2 matrix, each row of which 
            specifies a point in R^2. The 'Signature' field contains 
            the corresponding Euclidean signature. The 'Arcs' field 
            contains the bivertex arcs of the boundary as well as the 
            corresponding arcs of the Euclidean signature. The 
            'Pt2Arc' field contains a vector whose jth entry is the 
            bivertex arc to which the jth point of the boundary curve 
            belong, where a 0 indicates that the point does not 
            belong to a bivertex arc.
            
'placements':     This output is a 1-by-n structure whose jth 
            column corresponds to the jth piece placed by the 
            algorithm. The 'Piece' field indicates the index of the 
            piece within the 'pieces' structure. The 'Score' field 
            contains a matrix [q q_tilde] where these are the scores 
            resulting from the piece locking that placed the jth 
            piece as described in Sect. 5.3 of [1]. The 'g_lock' field
            contains the rigid motion that takes the piece boundary 
            to its place in the solution, represented as a matrix 
            [theta a b] where these parameters are as described in 
            Sect. 3.3 of [2]. The 'Fit' field contains a structure 
            containing detailed information about the fit that was 
            refined by piece-locking to produce 'g_lock'. The field
            'Neighbors' contains the indices (within the 'pieces'
            structure) of the pieces adjacent to the jth piece in the
            solution.
            
'tracker':  This ouput is a 1-by-n structure whose jth column gives
            detailed information about the state of the puzzle 
            assembly after j pieces have been placed. The fields 
            'PPc', 'RPc', 'IArcs', 'AArcs', 'IPts', 'APts', and 
            'Pc2Place' indicate, respectively, the placed pieces, the
            remaining pieces, the inactive arcs, the active arcs, the 
            inactive points, the active points, and order in which 
            pieces have been placed.

'parameters':   This variable contains the parameters (described in 
            [1]) used by the algorithm. Use the Solve_Puzzle() 
            function for a GUI to input these parameters, or see its 
            code for the form this input must take.

%--------------------------------------------------------------------
%}



function [pieces placements tracker parameters] = Assemble(puzzle, plotter, saver, parameters, placements, tracker)

%--------------------------------------------------------------------------
% Preliminaries
%--------------------------------------------------------------------------

% Check to see if parameters are appropriate
errors = check_parameters(parameters, 1);
if(~isempty(errors))
    msgbox(errors);
    return;
end;

% Isolate parameters
alpha = parameters{1};
beta = parameters{2};
gamma = parameters{3};
C_1 = parameters{4};
C_2 = parameters{5};
K_1 = parameters{6};
K_2 = parameters{7};
K_4 = parameters{8};
lambda_0 = parameters{9};
lambda_1 = parameters{10};
nu = parameters{11};
epsilon = parameters{12};
rho = parameters{13};
j_max = parameters{14};
eta_1 = parameters{15};
eta_2 = parameters{16};
Q_1 = parameters{17};
Q_2 = parameters{18};
Q_2_star = parameters{19};
Q_3 = parameters{20};
parameter_sequence = parameters{21};
j_star = size(parameter_sequence, 2);



% Set up figure
if(plotter > 1)
    fh = figure('Units', 'normalized', 'OuterPosition', [.225, .1, .45, .8], 'Name', 'Current Subpuzzle');
else
    fh = [];
end;

% Set up progress report
[handles ps] = progress_report({'Calculating Euclidean Signatures:'; 'Approximating Bivertex Arc Decompositions:' ; 'Placing Pieces:' ;...
    {'Comparing Pieces Arc by Arc:'  ; 'Checking Fits'}});

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Initialize variables as necessary
%--------------------------------------------------------------------------
if(~isa(puzzle, 'struct'))
    nPieces = size(puzzle, 2);
    
    % Create pieces structure
    pieces = struct('Points', {}, 'Signature', {}, 'Arcs', {}, 'Pt2Arc', {}, 'Active', {});
    for c1 = 1:nPieces    
        pieces(1, c1).Points = puzzle{1, c1};
    end;

    % Calculate the Euclidean signatures of pieces
    handles = update_report(handles, 1, 1);
    weights = zeros(1, nPieces);
    av_size = 0;
    av_len = 0;
    for c1 = 1:nPieces

        % Open curves if closed
        if(pieces(c1).Points(end, :) == pieces(c1).Points(1, :))
            pieces(c1).Points(end, :) = [];
        end;

        % Calculate Euclidean Signatures
        y = Euclidean_Signature(pieces(c1).Points);  

        % Orient curves
        if(sum(y(:, 1)) < 0)
            y = Orient_Reverse(y);
            pieces(c1).Points = flipud(pieces(c1).Points);
            disp(['Orientation of curve ' num2str(c1) ' reversed']);
        end;

        % Add signature to variable "Pieces"
        pieces(c1).Signature = y;

        % Calculate weight
        weights(1, c1) = sum(abs(y(:, 1)));

        % Update progress report
        ps{1, 1} = c1/nPieces;
        handles = update_report(handles, ps);

        av_size = av_size + size(pieces(c1).Points, 1);
        av_len = av_len + sum(vecnorm(pieces(c1).Points - circshift(pieces(c1).Points, [-1, 0])));
    end;
    av_size = av_size/nPieces;
    av_len = av_len/nPieces;
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------
    % Split pieces into arcs
    %--------------------------------------------------------------------------

    % Update progress report
    handles = update_report(handles, 2, 1);

    % Calculate characteristic distances
    [D_x D_y] = max_d({pieces.Points});
    [D_kappa D_kappa_s] = max_d({pieces.Signature});

    % Account for trivial cases
    if D_x == 0
        D_x = 1;
    end;   
    if D_y == 0
        D_y = 1;
    end;   
    if D_kappa == 0
        D_kappa = 1;
    end;   
    if D_kappa_s == 0
        D_kappa_s = 1;
    end;   

    % Calculate cutoffs
    d0 = D_kappa_s/lambda_0;
    d1 = D_kappa/lambda_1;

    % Calculate approximate bivertex arc decompositions
    for c1 = 1:nPieces
        [pieces(c1).Arcs pieces(c1).Pt2Arc] = Bivertex_Arc_Decomposition(pieces(c1).Points, pieces(c1).Signature, d0, d1);    

        %Update progress report
        ps{2, 1} = c1/nPieces;
        handles = update_report(handles, ps);
    end;
else
    pieces = puzzle;
    nPieces = size(puzzle, 2);
    
    % Calculate characteristic distances
    [D_x D_y] = max_d({pieces.Points});
    [D_kappa D_kappa_s] = max_d({pieces.Signature});

    % Account for trivial cases
    if D_x == 0
        D_x = 1;
    end;   
    if D_y == 0
        D_y = 1;
    end;   
    if D_kappa == 0
        D_kappa = 1;
    end;   
    if D_kappa_s == 0
        D_kappa_s = 1;
    end;  
    
    % Calculate averages
    av_size = 0;
    av_len = 0;
    for c1 = 1:nPieces
        av_size = av_size + size(pieces(c1).Points, 1);
        av_len = av_len + sum(vecnorm(pieces(c1).Points - circshift(pieces(c1).Points, [-1, 0])));
    end;
    av_size = av_size/nPieces;
    av_len = av_len/nPieces;
    
    
    % Update progress report
    handles = update_report(handles, 1, 1);
    ps{1, 1} = 1;
    handles = update_report(handles, ps);
    handles = update_report(handles, 2, 1);
    ps{2, 1} = 1;
    handles = update_report(handles, ps);
end;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Assemble puzzle piece by piece
%--------------------------------------------------------------------------

% Initialize variables if necessary
if(~exist('placements', 'var') || ~exist('tracker', 'var'))
% Choose starting piece
[dummy piece1] = max(weights);
pieces(piece1).Active = 1;
placements = struct('Piece', {piece1}, 'Score', {[]}, 'g_lock', {[0 0 0]}, 'Fit', {[]}, 'Neighbors', {[]});
tracker = struct('PPc', {piece1}, 'RPc', {1:nPieces}, 'IArcs', {cell(nPieces, 1)}, 'AArcs', {cell(nPieces, 1)},...
    'IPts', {cell(nPieces, 1)}, 'APts', {cell(nPieces, 1)}, 'Pc2Place', {zeros(nPieces, 1)}, 'SP_Bdry', {pieces(piece1).Points},...
    'SP_BdryPt_2_PcPt', {[ones(size(pieces(piece1).Points, 1), 1)*piece1 (1:size(pieces(piece1).Points, 1))']});
for c1 = 1:nPieces
    tracker(1).AArcs{c1, 1} = (1:size(pieces(c1).Arcs, 1))';
    tracker(1).APts{c1, 1} = (1:size(pieces(c1).Points, 1))';
end;
tracker(1).Pc2Place(piece1, 1) = 1;
tracker(1).RPc(piece1) = [];
end;

% Initialize plot if applicable
if(plotter > 1)
    figure(fh);
    farcpts1 = zeros(0, 2);
    farcpts2 = zeros(0, 2);
    ah = axes('ytick', [], 'xtick', []);
    title('Current Subpuzzle', 'FontUnits', 'normalized', 'FontSize', .04);
    axis equal;
    hold on;
    farch1 = [];
    farch2 = [];
    fragh = plot(tracker(end).SP_Bdry(:, 1), tracker(end).SP_Bdry(:, 2), 'xb'); 
    legend(fragh, 'Active Points', 'Location', 'Best');
    for c1 = 1:size(placements, 2)
        x2 = pieces(placements(c1).Piece).Points;
        x2 = transf(x2, placements(c1).g_lock);
        [txy] = mean(x2);
        text(txy(1, 1), txy(1, 2), num2str(c1));
        x2 = [x2 ; x2(1, :)];
        plot(x2(:, 1), x2(:, 2), 'Color', 'k', 'LineWidth', 1.1) 
    end;
end;

% Initialize variables
j = 1;
c1 = size(placements, 2)+1;
nRPc = nPieces - c1 + 2;
PScores = cell(nPieces, nPieces);
% Update progress report
handles = update_report(handles, 3, 1);
% Place pieces one by one
while(c1 < nPieces + 1)
    
    % Initialize Fits variable
    Fits = struct('Pieces', {}, 'Arcs', {}, 'ArcTrans', {}, 'g_fit', {}, 'Size', {}, 'Score', {}, 'Slot', {});
    
    % Update progress report
    handles = update_report(handles, [4 1], 1);
    
    % Examine all pairings of unplaced and placed pieces
    for c2 = 1:nPieces-c1+1
        for c3 = 1:c1-1
            %--------------------------------------------------------------------------
            % Find P-Scores between arcs of RePieces and arcs of AsPieces
            %--------------------------------------------------------------------------
            if(isempty(PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}))
                PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)} = zeros(size(pieces(tracker(end).RPc(1, c2)).Arcs, 1), size(pieces(tracker(end).PPc(1, c3)).Arcs, 1));
                for c4 = 1:size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1)
                    for c5 = 1:size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1)
                        PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)) = ...
                            Signature_Similarity_Coefficient(pieces(tracker(end).RPc(1, c2)).Arcs{tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1), 2},...
                            Orient_Reverse(pieces(tracker(end).PPc(1, c3)).Arcs{tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1), 2}), D_kappa, gamma, C_1, alpha, epsilon);
                    end;
                end;
                PScores{tracker(end).PPc(1, c3), tracker(end).RPc(1, c2)} = PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)};
            end;
            %--------------------------------------------------------------------------
           
            %--------------------------------------------------------------------------
            % Find sequences of consecutive high P-Scores
            %--------------------------------------------------------------------------
            included = zeros(size(PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}));
            for c4 = 1:size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1)
                for c5 = 1:size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1)
                    if(~included(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1))...
                            && PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)) >= parameter_sequence(j).p_0)
                        included(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)) = 1;
                        tArcs = [tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1) tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)];
                        c6 = 1;
                        while(PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}(tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4+c6-1, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5-c6-1, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1)) >= parameter_sequence(j).p_0)
                            if(abs(tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4+c6-1, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1)- tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4+c6-2, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1)) > 1 ...
                                    || abs(tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5-c6-1, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1) - tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5-c6-1+1, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1)) > 1)
                                break;
                            else
                                tArcs = [tArcs ; mod(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1)+c6-1, size(included, 1))+1  mod(tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)-c6-1, size(included, 2))+1];
                                included(mod(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1)+c6-1, size(included, 1))+1, mod(tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)-c6-1, size(included, 2))+1) = 1;
                                c6 = c6 + 1;
                            end;
                        end;
                        c7 = 1;
                        while(PScores{tracker(end).RPc(1, c2), tracker(end).PPc(1, c3)}(tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4-c7-1, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1), tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5+c7-1, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1)) >= parameter_sequence(j).p_0)                             
                            if(abs(tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4-c7-1, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1)- tracker(end).AArcs{tracker(end).RPc(1, c2)}(mod(c4-c7-1+1, size(tracker(end).AArcs{tracker(end).RPc(1, c2)}, 1))+1, 1)) > 1 ...
                                    || abs(tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5+c7-1, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1) - tracker(end).AArcs{tracker(end).PPc(1, c3)}(mod(c5+c7-2, size(tracker(end).AArcs{tracker(end).PPc(1, c3)}, 1))+1, 1)) > 1)
                                break;
                            else
                                tArcs = [mod(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1)-c7-1, size(included, 1))+1 mod(tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)+c7-1, size(included, 2))+1 ;tArcs];
                                included(mod(tracker(end).AArcs{tracker(end).RPc(1, c2)}(c4, 1)-c7-1, size(included, 1))+1, mod(tracker(end).AArcs{tracker(end).PPc(1, c3)}(c5, 1)+c7-1, size(included, 2))+1) = 1;
                                c7 = c7 + 1; 
                            end;
                        end;
                        if(c6+c7-1-min(size(included)) > 0)
                            for c8 = 0:c6+c7-1-min(size(included))
                                Fits(end+1).Size = min(size(included));
                                Fits(end).Pieces = [tracker(end).RPc(1, c2) tracker(end).PPc(1, c3)];
                                Fits(end).Arcs = tArcs(1+c8:Fits(end).Size+c8, :);
                                Fits(end).Slot = 1;
                            end;
                        else
                            Fits(end+1).Size = c6+c7-1;
                            Fits(end).Pieces = [tracker(end).RPc(1, c2) tracker(end).PPc(1, c3)];
                            Fits(end).Arcs = tArcs;
                            Fits(end).Slot = 1;
                        end;
                        if(Fits(end).Size < parameter_sequence(j).m_0)
                            Fits(end) = [];
                        end;
                    end;                    
                end;
            end;
            %--------------------------------------------------------------------------
            
            %Update progress report
            ps{4, 1}{1, 1} = ((c1-1)*(c2-1)+c3)/((c1-1)*(nPieces-c1+1));
            handles = update_report(handles, ps);
        end;
    end;
  
    %Sort fits by size
    [dummy FitIx] = sort([Fits.Size], 'descend');
    Fits = Fits(FitIx);
    
    % Initialize variables
    progress = 0;
    c2 = 1;
    % Update progress report
    handles = update_report(handles, [4 2], 1);
    % Check fits
    while(c2 < size(Fits, 2)+1)     
        if(pieces(Fits(c2).Pieces(1, 2)).Active && ismember(Fits(c2).Pieces(1, 1), tracker(end).RPc) && ~any(ismember(Fits(c2).Arcs(:, 2), tracker(end).IArcs{Fits(c2).Pieces(1, 2)})))
            
            %--------------------------------------------------------------------------
            % Calculate transformation and mu-score
            %--------------------------------------------------------------------------
            if(isempty(Fits(c2).Score))
                Fits(c2).ArcTrans = zeros(size(Fits(c2).Arcs, 1), 3);
                Fits(c2).g_fit = [0 0 0];
                for c3 = 1:Fits(c2).Size
                    Fits(c2).ArcTrans(c3, 1) = Rigid_Motion_Angle(pieces(Fits(c2).Pieces(1, 1)).Arcs{Fits(c2).Arcs(c3, 1), 1},...
                        flipud(transf(pieces(Fits(c2).Pieces(1, 2)).Arcs{Fits(c2).Arcs(c3, 2), 1}, placements(tracker(end).Pc2Place(Fits(c2).Pieces(1, 2))).g_lock)),...
                        pieces(Fits(c2).Pieces(1, 1)).Arcs{Fits(c2).Arcs(c3, 1), 2},...
                        Orient_Reverse(pieces(Fits(c2).Pieces(1, 2)).Arcs{Fits(c2).Arcs(c3, 2), 2}), beta);
                    Fits(c2).g_fit(1, 1) = Fits(c2).g_fit(1, 1) + cos(Fits(c2).ArcTrans(c3, 1)) + sqrt(-1)*sin(Fits(c2).ArcTrans(c3, 1));
                end;
                Fits(c2).g_fit(1, 1) = angle(Fits(c2).g_fit(1, 1));
                
                %Calculate translations with after rotation through average angle theta
                for c3 = 1:Fits(c2).Size  
                    [Fits(c2).ArcTrans(c3, 2) Fits(c2).ArcTrans(c3, 3)] = Rigid_Motion_Translation(pieces(Fits(c2).Pieces(1, 1)).Arcs{Fits(c2).Arcs(c3, 1), 1},...
                        flipud(transf(pieces(Fits(c2).Pieces(1, 2)).Arcs{Fits(c2).Arcs(c3, 2), 1}, placements(tracker(end).Pc2Place(Fits(c2).Pieces(1, 2))).g_lock)),...
                        pieces(Fits(c2).Pieces(1, 1)).Arcs{Fits(c2).Arcs(c3, 1), 2},...
                        Orient_Reverse(pieces(Fits(c2).Pieces(1, 2)).Arcs{Fits(c2).Arcs(c3, 2), 2}), Fits(c2).g_fit(1, 1), beta);
                end;   

                Fits(c2).Score = muScore(Fits(c2).ArcTrans, D_x, D_y, C_2);
                Fits(c2).g_fit(1, 2:3) = mean(Fits(c2).ArcTrans(:, 2:3), 1);
            end;
            

            if(Fits(c2).Score < parameter_sequence(j).mu_0)
                % In this case, it may be that a subset of this collection
                % of arc pairings will have a high enough mu-score.        
                % We need only check this possibility if trimming an arc
                % will not drop the number of arc pairs below m_0
                if(Fits(c2).Size > parameter_sequence(j).m_0)
                    
                    % Create a temporary structure
                    tFits = struct('Pieces', {}, 'Arcs', {}, 'ArcTrans', {}, 'g_fit', {}, 'Size', {}, 'Score', {}, 'Slot', {});

                    % We generate smaller sequences of matched pairings by
                    % deleting the first and last arcs of Fits(c2). Note,
                    % however, that we only consider the result of deleting
                    % the last arc if Fits(c2) did not arise from deleting
                    % the first arc of a larger fit (i.e., if Fits(c2).Slot
                    % == 1). This avoids generating duplicate fits. 
                    
                    % Delete last arc
                    if(Fits(c2).Slot)                
                        tFits(1).Pieces = Fits(c2).Pieces;               
                        tFits(1).Arcs = Fits(c2).Arcs(1:end-1, :);
                        tFits(1).Size = Fits(c2).Size - 1;
                        tFits(1).ArcTrans = Fits(c2).ArcTrans(1:end-1, :);               
                        tang = 0;
                        for c4 = 1:tFits(1).Size
                            tang = tang + cos(tFits(1).ArcTrans(c4, 1)) + sqrt(-1)*sin(tFits(1).ArcTrans(c4, 1));
                        end;
                        tang = angle(tang);
                        tFits(1).g_fit = [tang 0 0];
                        for c4 = 1:tFits(1).Size
                            [tFits(1).ArcTrans(c4, 2) tFits(1).ArcTrans(c4, 3)] = Rigid_Motion_Translation(pieces(tFits(1).Pieces(1, 1)).Arcs{tFits(1).Arcs(c4, 1), 1},...
                            flipud(transf(pieces(tFits(1).Pieces(1, 2)).Arcs{tFits(1).Arcs(c4, 2), 1}, placements(tracker(end).Pc2Place(tFits(1).Pieces(1, 2))).g_lock)),...
                            pieces(tFits(1).Pieces(1, 1)).Arcs{tFits(1).Arcs(c4, 1), 2},...
                            Orient_Reverse(pieces(tFits(1).Pieces(1, 2)).Arcs{tFits(1).Arcs(c4, 2), 2}), tFits(1).g_fit(1, 1), beta);
                        end;
                        tFits(1).Score = muScore(tFits(1).ArcTrans, D_x, D_y, C_2);
                        tFits(1).g_fit(1, 2:3) = mean(tFits(1).ArcTrans(:, 2:3), 1);
                        tFits(1).Slot = 1;
                    end;

                    % Delete first arc
                    tFits(end+1).Pieces = Fits(c2).Pieces;
                    tFits(end).Arcs = Fits(c2).Arcs(2:end, :);
                    tFits(end).Size = Fits(c2).Size - 1;
                    tFits(end).ArcTrans = Fits(c2).ArcTrans(2:end, :);
                    tang = 0;
                    for c4 = 1:tFits(end).Size
                        tang = tang + cos(tFits(end).ArcTrans(c4, 1)) + sqrt(-1)*sin(tFits(end).ArcTrans(c4, 1));
                    end;
                    tang = angle(tang);               
                    tFits(end).g_fit = [tang 0 0];
                    for c4 = 1:tFits(end).Size
                        [tFits(end).ArcTrans(c4, 2) tFits(end).ArcTrans(c4, 3)] = Rigid_Motion_Translation(pieces(tFits(end).Pieces(1, 1)).Arcs{tFits(end).Arcs(c4, 1), 1},...
                        flipud(transf(pieces(tFits(end).Pieces(1, 2)).Arcs{tFits(end).Arcs(c4, 2), 1}, placements(tracker(end).Pc2Place(tFits(end).Pieces(1, 2))).g_lock)),...
                        pieces(tFits(end).Pieces(1, 1)).Arcs{tFits(end).Arcs(c4, 1), 2},...
                        Orient_Reverse(pieces(tFits(end).Pieces(1, 2)).Arcs{tFits(end).Arcs(c4, 2), 2}), tFits(end).g_fit(1, 1), beta);
                    end;
                    tFits(end).Score = muScore(tFits(end).ArcTrans, D_x, D_y, C_2);
                    tFits(end).g_fit(1, 2:3) = mean(tFits(end).ArcTrans(:, 2:3), 1);
                    tFits(end).Slot = 0;
              
                    
                    % The new fits in tFits are shorter, so we must figure 
                    % out where to place them within Fits in order to 
                    % preserve the sorting on that structure 
                    c3 = c2 + 1;
                    while(c3 < size(Fits, 2) + 1 && Fits(c3).Size == Fits(c2).Size)                   
                        c3 = c3 + 1;
                    end;
                    
                    % Insert tFits into Fits. 
                    if(c3 <= size(Fits, 2))
                        Fits = [Fits(1:c3-1) tFits Fits(c3:end)];
                    else
                        Fits = [Fits(1:c3-1) tFits];
                    end;
                
                end;
                % Delete the fit that has now been subdivided.
                Fits(c2) = []; 
                % Update counter variable
                c2 = c2 - 1;
                
            else
                % In this case we perform piece locking to refine the
                % transformation and provide additional confidence scores
                           
                [g_lock tPiecePtIcs_3 tracker(end).SP_Bdry_PtIcs_3 tPiecePtIcs tracker(end).SP_Bdry_PtIcs xph] =...
                    Lock(Fits(c2).g_fit, pieces(Fits(c2).Pieces(1, 1)).Points, tracker(end).SP_Bdry, K_1, K_2, parameter_sequence(j).K_3, K_4, epsilon, nu, rho, j_max, (plotter > 2), fh);
                
				% Compute Scores
				q_1 = size(tracker(end).SP_Bdry_PtIcs_3, 1)/size(tracker(end).SP_Bdry_PtIcs, 1);
                q_2 = 0;
                tn = size(tPiecePtIcs_3, 1);
                for c3 = 1:tn
                    if(abs(tPiecePtIcs_3(mod(c3+1-1, tn)+1)-tPiecePtIcs_3(c3)) == 1 || abs(tPiecePtIcs_3(mod(c3+1-1, tn)+1)-tPiecePtIcs_3(c3)) == size(pieces(Fits(c2).Pieces(1, 1)).Points, 1)-1)
                        q_2 = q_2 + vecnorm(pieces(Fits(c2).Pieces(1, 1)).Points(tPiecePtIcs_3(mod(c3+1-1, tn)+1), :)-pieces(Fits(c2).Pieces(1, 1)).Points(tPiecePtIcs_3(c3), :));
                    end;
                end;
                q_2 = q_2/av_len;
				q_3 = sum(abs(pieces(Fits(c2).Pieces(1, 1)).Signature(tPiecePtIcs_3, 1)))/(D_kappa*av_size);
				
                % Update progress report
                ps{4, 1}{2, 1} = c2/(size(Fits, 2)+1);
                handles = update_report(handles, ps);
                
                % If scores are high enough, place the piece
                if((eta_1*q_1+eta_2*q_3 > Q_1 || q_2 > Q_2_star) && q_2 > Q_2 && q_3 > Q_3)                   
                    
                    % Add new placement to placements variable
                    placements(c1).Score = [q_1 q_2 q_3 eta_1*q_1+eta_2*q_3];
                    placements(c1).Piece = Fits(c2).Pieces(1, 1);
                    placements(c1).g_lock = g_lock;
                    placements(c1).Fit = Fits(c2);
                    
                    % Update tracker variable
                    tracker(end+1) = tracker(end);
                    tracker(end).Pc2Place(placements(c1).Piece, 1) = c1;                   
                    
    
                    %Mark arcs/points as used or remaining and introduce neighbors/activate
                    for c3 = 1:size(tPiecePtIcs, 1)                      
                        if(pieces(Fits(c2).Pieces(1, 1)).Pt2Arc(tPiecePtIcs(c3, 1), 1))
                            tracker(end).IArcs{Fits(c2).Pieces(1, 1), 1}(end+1, 1) = pieces(Fits(c2).Pieces(1, 1)).Pt2Arc(tPiecePtIcs(c3, 1), 1);
                        end;                     
                    end;
                    placements(c1).Neighbors = [];
					pieces(Fits(c2).Pieces(1, 1)).Active = 2;
                    for c3 = 1:size(tracker(end).SP_Bdry_PtIcs, 1)
                        % Introduce new neighbors
                        if(~ismember(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1), placements(c1).Neighbors))
                            placements(c1).Neighbors = [placements(c1).Neighbors tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1)];
                            pieces(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1)).Active = 2;
							[pieces(placements(tracker(c1).Pc2Place(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1))).Neighbors).Active] = deal(2);
							placements(tracker(c1).Pc2Place(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1))).Neighbors = [placements(tracker(c1).Pc2Place(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1))).Neighbors Fits(c2).Pieces(1, 1)];
                        end;
                        
                        % Mark arcs as inactive
                        if(pieces(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1)).Pt2Arc(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 2), 1))
                            tracker(end).IArcs{tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1), 1}(end+1, 1) = pieces(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1)).Pt2Arc(tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 2), 1);
                        end;
                        
                        %Mark points as inactive
                        tracker(end).IPts{tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1), 1} = [tracker(end).IPts{tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 1), 1} ; tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs(c3, 1), 2)];
                    end;
                    tracker(end).IPts{Fits(c2).Pieces(1, 1), 1} = tPiecePtIcs;
                    tracker(end).APts{Fits(c2).Pieces(1, 1), 1} = (1:size(pieces(Fits(c2).Pieces(1, 1)).Points, 1))';
                    tracker(end).APts{Fits(c2).Pieces(1, 1), 1}(tracker(end).IPts{Fits(c2).Pieces(1, 1), 1}, :) = [];                   
                    tracker(end).AArcs{Fits(c2).Pieces(1, 1), 1} = (1:size(pieces(Fits(c2).Pieces(1, 1)).Arcs, 1))'; 
                    tracker(end).AArcs{Fits(c2).Pieces(1, 1), 1}(tracker(end).IArcs{Fits(c2).Pieces(1, 1), 1}, :) = [];                   
                    for c3 = 1:c1-1
                        %Calculate AArcs
                        tracker(end).AArcs{tracker(end).PPc(1, c3), 1} = (1:size(pieces(tracker(end).PPc(1, c3)).Arcs, 1))'; 
                        tracker(end).AArcs{tracker(end).PPc(1, c3), 1}(tracker(end).IArcs{tracker(end).PPc(1, c3), 1}, :) = [];
                        
                        %Calculate APts
                        tracker(end).APts{tracker(end).PPc(1, c3), 1} = (1:size(pieces(tracker(end).PPc(1, c3)).Points, 1))';
                        tracker(end).APts{tracker(end).PPc(1, c3), 1}(tracker(end).IPts{tracker(end).PPc(1, c3), 1}, :) = [];
                    end;

                    
                    %Calculate new Fragment
                    tracker(end).SP_Bdry(tracker(end).SP_Bdry_PtIcs, :) = [];
                    tPiece = pieces(Fits(c2).Pieces(1, 1)).Points;
                    tTrack = (1:size(tPiece, 1))';                
                    tPiece(tPiecePtIcs, :) = [];
                    tTrack(tPiecePtIcs, :) = [];
                    tracker(end).SP_Bdry = [tracker(end).SP_Bdry ; transf(tPiece, g_lock)];
                    tracker(end).SP_BdryPt_2_PcPt(tracker(end).SP_Bdry_PtIcs, :) = [];
                    tracker(end).SP_BdryPt_2_PcPt = [tracker(end).SP_BdryPt_2_PcPt ; [ones(size(tTrack, 1), 1)*placements(c1).Piece tTrack]];
                    
                    %Mark pieces as used or remaining
                    tracker(end).PPc(1, end+1) = placements(c1).Piece;
                    tracker(end).RPc = 1:nPieces;
                    tracker(end).RPc(tracker(end).PPc) = [];
                      
                    %Update progress report
                    ps{3, 1} = (c1-nPieces+nRPc)/nRPc;
                    handles = update_report(handles, ps);

                    %Display, if appropriate
                    if(plotter > 1)
                        figure(fh);
                        axis auto;
                        
                        %Update Fragment
                        delete(fragh);
                        fragh = plot(tracker(end).SP_Bdry(:, 1), tracker(end).SP_Bdry(:, 2), 'xb');
                        
                        %Plot new piece
                        x2 = pieces(placements(c1).Piece).Points;
                        x2 = transf(x2, placements(c1).g_lock);
                        [txy] = mean(x2);
                        text(txy(1, 1), txy(1, 2), num2str(c1));
                        x2 = [x2 ; x2(1, :)];
                        plot(x2(:, 1), x2(:, 2), 'Color', 'k', 'LineWidth', 1.1)
                        
                        %Update fit arc points
                        for c4 = 1:size(placements(c1).Fit.Arcs, 1) 
                            x1 = pieces(placements(c1).Fit.Pieces(1, 1)).Arcs{placements(c1).Fit.Arcs(c4, 1), 1};
                            x1 = transf(x1, placements(tracker(end).Pc2Place(placements(c1).Fit.Pieces(1, 1), 1)).g_lock);
                            x2 = pieces(placements(c1).Fit.Pieces(1, 2)).Arcs{placements(c1).Fit.Arcs(c4, 2), 1};
                            x2 = transf(x2, placements(tracker(end).Pc2Place(placements(c1).Fit.Pieces(1, 2), 1)).g_lock);
                            farcpts1 = [farcpts1 ; x1];
                            farcpts2 = [farcpts2 ; x2];
                        end;
                        delete(farch1);
                        delete(farch2);
                        farch1 = plot(farcpts1(:, 1), farcpts1(:, 2), 'm*', 'MarkerSize', 3);
                        farch2 = plot(farcpts2(:, 1), farcpts2(:, 2), 'r*', 'MarkerSize', 3);
                        legend([fragh farch1 farch2], 'Active Points', 'Matched Bivertex Arcs', 'Matched Bivertex Arcs', 'Location', 'Best');
                    end;
                    
                    % Indicate that progress was made 
                    progress = 1;
                    
                    % Save current data if applicable
                    if(saver)
                        save('CurrentPuzzleData', 'pieces', 'placements', 'tracker', 'parameters');
                    end;
                    
                    % Update counting variable
                    c1 = c1 + 1;
                end;
                
                % Delete the extra plot possibly generated by Lock(). We
                % wait until now in order to provide a more continuous
                % appearance when plotter = 3.
                delete(xph);
            end;

        end;
        
        % Update counting variable             
        c2 = c2 + 1;
    end;    
    % Update progress report
    handles = update_report(handles, [4 2], 0); 
    ps{4, 1}{2, 1} = 0;
    handles = update_report(handles, ps);
    
    
    % If the algorithm has dead-ended, increase depth in paraemter sequence
    % or terminate algorithm if this is not possible
    if(~progress)
        if(j < j_star)
            j = j + 1;
            %Activate pieces
            for c4 = 1:size(tracker(end).PPc, 2)
                pieces(tracker(end).PPc(1, c4)).Active = 1;
            end;
        else            
            display('Dead-ended Early');
            if(plotter == 1)    
                fh = figure('Units', 'normalized', 'OuterPosition', [.225, .1, .45, .8], 'Name', 'Puzzle');
                axes('ytick', [], 'xtick', []);
                hold on;
                axis equal;
                for c2 = 1:size(placements, 2)
                    x = pieces(placements(c2).Piece).Points;
                    x = transf(x, placements(c2).g_lock);
                    x = [x ; x(1, :)];
                    plot(x(:, 1), x(:, 2), 'Color', 'k', 'LineWidth', 1.1)
                end;  
            end;            
            return;
        end;      
    else
        j = 1;
        %Deactivate pieces
        for c4 = 1:size(tracker(end).PPc, 2)
            pieces(tracker(end).PPc(1, c4)).Active = max(0, pieces(tracker(end).PPc(1, c4)).Active-1);
        end;
    end;	
end;

% Plot final solution, if applicable
if(plotter == 1)    
    fh = figure('Units', 'normalized', 'OuterPosition', [.225, .1, .45, .8], 'Name', 'Puzzle');
    axes('ytick', [], 'xtick', []);
    hold on;
    axis equal;
    for c2 = 1:size(placements, 2)
        x = pieces(placements(c2).Piece).Points;
        x = transf(x, placements(c2).g_lock);
        x = [x ; x(1, :)];
        plot(x(:, 1), x(:, 2), 'Color', 'k', 'LineWidth', 1.1)
    end;  
end;
%--------------------------------------------------------------------------


end  


%--------------------------------------------------------------------
%--------------------------------------------------------------------
%--------------------------------------------------------------------
% Auxiliary Functions
%--------------------------------------------------------------------
%--------------------------------------------------------------------
%--------------------------------------------------------------------



function new_points = transf(points, trans)

% This function applies a rigid motion to inputted points
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'points':   This should be n points represented as an n-by-2 matrix, 
            each row of which specifies a point in R^2. 

'trans':    This input should be a matrix of the form [theta a b] 
            where these parameters specify a rigid motion as 
            described in [2] (a translation by [a b] and a rotation 
            by theta radians around the origin). 

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'new_points':   This output gives the transformed points as an n-by-2 
            matrix, each row of which specifies a point in R^2.

%--------------------------------------------------------------------
%}

new_points(:, 1) = cos(trans(1))*points(:, 1) - sin(trans(1))*points(:, 2)+trans(2);
new_points(:, 2) = sin(trans(1))*points(:, 1) + cos(trans(1))*points(:, 2)+trans(3);


end


%--------------------------------------------------------------------
%--------------------------------------------------------------------

    
function [hspread vspread] = max_d(curves)

% This function computes the greatest horizontal and vertical distances
% present in a collection of curves
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'points':   This should be a single row cell array each element of 
            which contains a curve represented as an n-by-2 matrix, 
            each row specifying a point in R^2. 

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'hspread':  This output gives the largest horizontal distance present
            in any of the members of 'curves.'

'vspread':  This output gives the largest vertical distance present 
            in any of the members of 'curves.'

%--------------------------------------------------------------------
%}

    hspread = 0;
    vspread = 0;
    for c1 = 1:size(curves, 2)
        x = curves{1, c1};
        hspread = max(max(x(:, 1))-min(x(:, 1)), hspread);
        vspread = max(max(x(:, 2))-min(x(:, 2)), vspread);
    end;
                        
end  


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function new_signature = Orient_Reverse(signature)

% This function changes a Euclidean signature to reflect a reverse of
% orientation in the parent curve
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'signature':    This should be a discrete approximation to a 
            Euclidean signature, represented as an n-by-2 matrix, 
            each row of which specifies a point in R^2. The curve is 
            assumed to be closed; no repetition of points is 
            necessary. For example, 'signature' could be the output
            of Euclidean_Signature() called on a curve.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'new_signature':    This output gives the discrete Euclidean 
            signature for the curve optained by reversing the 
            orientation of the curve whose Euclidean signature is 
            given in 'signature.'

%--------------------------------------------------------------------
%}

signature(:, 1) = -signature(:, 1) ;
signature = flipud(signature);
new_signature = signature;

end


%--------------------------------------------------------------------
%--------------------------------------------------------------------



function [BA Pt2Arc] = Bivertex_Arc_Decomposition(curve, signature, d_0, d_1)


% This function approximates a bivertex decomposition as described in 
% Sect. 3.1 of [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'curve':    This should be a discretized planar curve of n points 
            represented as an n-by-2 matrix, each row of which 
            specifies a point in R^2. The curve is assumed to be 
            closed; no repetition of points is necessary.

'signature':This should be a discrete approximation to the Euclidean
            signature of 'curve' represented as an n-by-2 matrix, 
            each row of which specifies a point in R^2. The curve is 
            assumed to be closed; no repetition of points is 
            necessary. For example, 'signature' could be the output
            of Euclidean_Signature() called on 'curve'.

'd_0':      This positive real number specifies a kappa_s cut-off as
            described in Sect 3.1.

'd_1':      This postive real number specifies a minimum meaningful 
            change in kappa as described in Sect 3.1

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'BA':       This output gives an approximate bivertex decomposition 
            for 'curve' represented as an nA-by-2 cell array, where
            nA is the number of arcs in the decomposition. Each row
            of 'BA' corresponds to an arc; the first column contains 
            correponding points from 'curve' and the second column
            contains corresponding points from 'signature', 
            repsented as two-column matrices (as inputted).

'Pt2Arc':   This output is an n-by-1 matrix whose jth row contains 
            the index (within 'BA') of the bivertex arc to which the 
            jth point of 'curve' belongs, provided that the jth point 
            belongs to such an arc. If it does not, the value of jth 
            row is 0. 

%--------------------------------------------------------------------
%}


% Initialize variables
n1 = size(curve, 1);    
BA_indices = zeros(floor(n1/2), 2);


% Determine bivertex arc endpoint indices based on d_0 cut-off
if(abs(signature(1, 2)) < d_0 || signature(1, 2) == d_0)
    shifter = 0;
else
    shifter = 1;
end;
flipper = 0;
nA = 1;
for c1 = 1:n1
    if( (abs(signature(c1, 2))-d_0)*(abs(signature(mod(c1+1-1, n1)+1, 2))-d_0) <= 0)
        BA_indices(nA, mod(nA+shifter+1+flipper, 2) + 1) = c1;
        toggler = mod(nA+shifter+flipper+1,2);
        nA = nA + toggler;
        flipper = flipper+mod(toggler+1,2);
    else 
    if(abs(signature(c1,2)) > d_0 && abs(signature(mod(c1+1-1, n1)+1, 2)) > d_0 && (signature(c1, 2)*signature(mod(c1+1-1, n1)+1, 2)<0 ...
	|| (abs(signature(c1, 2))-abs(signature(mod(c1-1-1,n1)+1, 2)) < 0 && abs(signature(mod(c1+1-1, n1)+1, 2))-abs(signature(c1, 2)) > 0)) )
        BA_indices(nA, mod(nA+shifter+1+flipper, 2) + 1) = c1;
        toggler = mod(nA+shifter+flipper+1,2);
        nA = nA + toggler;
        flipper = flipper+mod(toggler+1,2);
        BA_indices(nA, mod(nA+shifter+1+flipper, 2) + 1) = c1;
        toggler = mod(nA+shifter+flipper+1,2);
        nA = nA + toggler;
        flipper = flipper+mod(toggler+1,2);
    end;
	end;
end;
if(shifter == 1)
    BA_indices(1, 1) = BA_indices(nA, 1);
end;

% Trim repeated arc
BA_indices(nA:end, :) = [];
nA = nA - 1;


% Reject arcs not meeting d_1 cut-off
c1 = 1;
while c1 < nA + 1           
    if(max(abs(signature(BA_indices(c1, 2), 1)), abs(signature(BA_indices(c1, 1), 1))) < d_1)              
        BA_indices(c1, :) = [];                
        nA = nA-1;
        c1 = c1 - 1;
    end;
    c1 = c1 + 1;
end;

% Assign ouput
BA = cell(nA,2);
Pt2Arc = zeros(n1, 1);
if(BA_indices(1, 1) > BA_indices(1, 2))
    BA{1, 1} = [curve(BA_indices(1, 1):end, :) ; curve(1: BA_indices(1, 2), :)];
    BA{1, 2} = [signature(BA_indices(1, 1):end, :) ; signature(1: BA_indices(1, 2), :)];
    Pt2Arc(BA_indices(1, 1) : end, 1) = 1;
    Pt2Arc(1:BA_indices(1, 2), 1) = 1;
else    
    BA{1, 1} = curve(BA_indices(1, 1):BA_indices(1, 2), :);
    BA{1, 2} = signature(BA_indices(1, 1):BA_indices(1, 2), :);
    Pt2Arc(BA_indices(1, 1) : BA_indices(1, 2), 1) = 1;
end;
for c1 = 2:nA          
    BA{c1, 1} = curve(BA_indices(c1, 1) : BA_indices(c1, 2), :);
    BA{c1, 2} = signature(BA_indices(c1, 1) : BA_indices(c1, 2), :);
    Pt2Arc(BA_indices(c1, 1) : BA_indices(c1, 2), 1) = c1;
end;

end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function p = Signature_Similarity_Coefficient(signature1, signature2, D, gamma, C_1, alpha, epsilon)

% This function calculates a signature similarity coefficient as 
% described in Sect. 3.2 of [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'signature1':   This should be a discrete approximation to a 
            Euclidean signature represented as an n-by-2 matrix, each 
            row of which specifies a point in R^2. For example, 
            'signature1' could be the output of the 
            Euclidean_Signature() function.

'signature2':   This should be a discrete approximation to a 
            Euclidean signature represented as an n-by-2 matrix, each 
            row of which specifies a point in R^2. For example, 
            'signature2' could be the output of the 
            Euclidean_Signature() function.

'D':        This positive real number specifies the scale of the
            comparison.

'gamma':    This real number specifies the exponent controlling the 
            simulated electro-static attraction.

'C_1':      This positive real number determines the distribution
            of p scores along the interval [0, 1].

'alpha':    This real number specifies the power by which curvature
            weights an average.

'epsilon':  This small postive real number serves as a cut-off to
            avoid infinities.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'p':        The signature similarity coeffient resulting from input.

%--------------------------------------------------------------------
%}


% Initialize variables
n1 = size(signature1, 1);
n2 = size(signature2, 1);    
h = zeros(n1, 1);

% Calculate strengths of correspondence
for i = 1:n1
    % Calculate and sum h over j
    for j = 1:n2
        d = norm(signature1(i, :) - signature2(j, :));
        if d < D 
            r = d/(D - d);
            h(i) = h(i) + 1/(r^gamma+epsilon);   
        end;
    end;
end;

% Rescale h for each i
p_point = h./(h + C_1);

% Take a weighted average over i
p_hat_1 = sum(p_point.*(abs(signature1(:, 1)).^alpha))/sum(abs(signature1(:, 1)).^alpha);

% Reverse roles and recompute
y3 = signature1;
signature1 = signature2;
signature2 = y3;
n3 = n1;
n1 = n2;
n2 = n3;
h = zeros(n1, 1);        

% Calculate strengths of correspondence
for i = 1:n1
    % Calculate and sum h over j
    for j = 1:n2
        d = norm(signature1(i, :) - signature2(j, :));
        if d < D  
            r = d/(D - d);
            h(i) = h(i) + 1/(r^gamma+epsilon);    
        end;
    end;
end;

% Rescale h for each i
p_point = h./(h + C_1);

% Take a weighted average over i
p_hat_2 = sum(p_point.*(abs(signature1(:, 1)).^alpha))/sum(abs(signature1(:, 1)).^alpha);


% Take minimum and assign output
p = min([p_hat_1, p_hat_2]);

 
end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function theta = Rigid_Motion_Angle(curve1, curve2, signature1, signature2, beta)

% This function reconstructs a rigid motion angle as described in 
% Sect. 3.3 of [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'curve1':   This should be a discretized planar curve of n points 
            represented as an n-by-2 matrix, each row of which 
            specifies a point in R^2. 

'curve2':   This should be a discretized planar curve of m points 
            represented as an m-by-2 matrix, each row of which 
            specifies a point in R^2. 

'signature1':   This should be a discrete approximation to the 
            Euclidean signature of 'curve1', represented as an 
            n-by-2 matrix, each row of which specifies a point in 
            R^2. 

'signature2':   This should be a discrete approximation to the 
            Euclidean signature of 'curve2', represented as an 
            m-by-2 matrix, each row of which specifies a point in 
            R^2. 

'beta':     This real number specifies the power by which kappa_s
            weights an average as described in Sect. 3.3.

%--------------------------------------------------------------------



%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'theta':    The resulting rigid motion angle, as described in
            Sect. 3.3. of [2].

%--------------------------------------------------------------------
%}



% Approximate tangent angles

tan1 = circshift(curve1, -1)-circshift(curve1, 1);
tan1(1, :) = tan1(2, :);
tan1(end, :) = tan1(end-1, :);
imtan1 = tan1(:, 1) + 1i*tan1(:,2);

tan2 = circshift(curve2, -1)-circshift(curve2, 1);
tan2(1, :) = tan2(2, :);
tan2(end, :) = tan2(end-1, :);
imtan2 = tan2(:, 1) + 1i*tan2(:,2);


% Find average tangent angle weighted by |\kappa_s|
angle1 = angle(sum(imtan1.*abs(signature1(:, 2)).^beta)/sum(abs(signature1(:, 2)).^beta));
angle2 = angle(sum(imtan2.*abs(signature2(:, 2)).^beta)/sum(abs(signature2(:, 2)).^beta));
theta = mod(angle2 - angle1, 2*pi);


end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function [a b] = Rigid_Motion_Translation(curve1, curve2, signature1, signature2, theta, beta)

% This function reconstructs a rigid motion translation as described 
% in Sect. 3.3 of [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'curve1':   This should be a discretized planar curve of n points 
            represented as an n-by-2 matrix, each row of which 
            specifies a point in R^2. 

'curve2':   This should be a discretized planar curve of m points 
            represented as an m-by-2 matrix, each row of which 
            specifies a point in R^2. 

'signature1':   This should be a discrete approximation to the 
            Euclidean signature of 'curve1', represented as an 
            n-by-2 matrix, each row of which specifies a point in 
            R^2. 

'signature2':   This should be a discrete approximation to the 
            Euclidean signature of 'curve2', represented as an 
            m-by-2 matrix, each row of which specifies a point in 
            R^2. 

'theta':    This should be the angle of the rigid motion as 
            outputted by Rigid_Motion_Angle().

'beta':     This real number specifies the power by which kappa_s
            weights an average as described in Sect. 3.3.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'a':        The resulting rigid motion horizontal translation as 
            described in Sect. 3.3 of [2].

'b':        The resulting rigid motion vertical translation as 
            described in Sect. 3.3 of [2].

%--------------------------------------------------------------------
%}

% Calculate center of mass weighted by kappa_s
cm1 = sum(curve1.*abs([signature1(:, 2) signature1(:, 2)]).^beta, 1)./sum(abs([signature1(:, 2) signature1(:, 2)]).^beta, 1);
cm2 = sum(curve2.*abs([signature2(:, 2) signature2(:, 2)]).^beta, 1)./sum(abs([signature2(:, 2) signature2(:, 2)]).^beta, 1);

% Calculate translation parameters
a = cm2(1, 1)-cm1(1, 1)*cos(theta)+cm1(1, 2)*sin(theta);
b = cm2(1, 2)-cm1(1, 1)*sin(theta)-cm1(1, 2)*cos(theta);

end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function mu = muScore(trans, D_x, D_y, C_2)

% This function calculates the mu score of a collection of rigid
% motions as described in [2].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'trans':    This should be a n-by-3 matrix each row of which 
            represents a transformation as [theta a b] where these 
            parameters specify a rigid motion as described in [2].

'D_x':      The characteristic horizontal distance used in the mu 
            score calculation as described in [2].

'D_y':      The characteristic vertical distance used in the mu score
            calculation as described in [2].

'C_2':      The parameter used in the mu score calculation as 
            described in [2].

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'mu':       This output gives resulting mu score as described in [2].

%--------------------------------------------------------------------
%}

% Calculate differences between angles
difang = zeros(size(trans, 1), size(trans, 1));
for c1 = 1:size(trans, 1) 
    for c2 = c1 + 1:size(trans, 1)
        difang(c1, c2) = mod(abs(trans(c1, 1) - trans(c2, 1)), 2*pi);
        if(difang(c1, c2) > pi)
            difang(c1, c2) = 2*pi - difang(c1, c2);
        end;
    end;
end;

% Calculate mu score
mu = 1-C_2*sum([max(max(difang))/pi  range(trans(:, 2:3))./[D_x  D_y]], 2);
if mu < 0
    mu = 0;
end;


end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function [handles ps] = progress_report(titles, namer)

% This function generates a progress report.
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'titles':   This variable should be a single column cell array each 
            of whose elements is a string or a further cell array of 
            the same form. A cell array element creates progress 
            bars whose progress determines the progress of the 
            immediately preceding bar. For instance, the input 
            "{'Calculating Euclidean Signatures:'; 'Approximating 
            Bivertex Arc Decompositions:' ; 'Placing Pieces:' ;
            {'Comparing Pieces Arc by Arc:'  ; 'Checking Fits'}}" 
            would generate a progress report with 5 bars, the 
            'Comparing Pieces Arc by Arc:' and 'Checking Fits' bars 
            showing progress within 'Placing Pieces'.

'namer':    This string will be the Name of the figure window that 
            contains the progress report. 

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'handles':  This output is a cell array containing the handles for 
            objects created in the progress report.

'ps':       This output is a cell array of the same form as 'titles', 
            but with elements of the interval [0, 1] in place of 
            strings, each corresponding to the progress of the 
            matching string. Since this function initializes the 
            progress report, these values all begin at 0.

%--------------------------------------------------------------------
%}

% Set default parameter value
if(~exist('namer', 'var'))
    namer = 'Progress Report';
end;

% Initialize variables
n = size(titles, 1);
fh = figure('Name', namer, 'Units', 'Points');
fdim = get(fh, 'Position');
w = fdim(1, 3);
timerhs = cell(n, 1);
ahs = cell(n, 1);
ths = cell(n, 2);
shs = cell(n, 1);
ps = cell(n, 1);
counter = 0;

% Create bars
for c1 = 0:n-1
    depth = 0;
    [timerhs{n-c1, 1} ahs{n-c1, 1} ths{n-c1, 1} shs{n-c1, 1} ps{n-c1, 1} w counter depth] = makebar(titles{n-c1, 1}, w, counter, depth);
end;

% Resize figure
set(fh, 'Position', [fdim(1, 1) fdim(1, 2) w counter*70+40]);
set(get(fh, 'Children'), 'Units', 'Normalized');
set(fh, 'Units', 'Normalized');
set(fh, 'OuterPosition', [.675 .1 .3 .8]);

% Set output
handles = {fh ahs ths shs timerhs};


end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function [timerhs ahs ths shs ps w counter depth] = makebar(names, w, counter, depth)  

% This recursively defined function is for use by the progress_report()
% function.

    if(isa(names, 'char'))
        counter = counter + 1;
        ps = 0;
        timerhs = [];
        ahs = axes('Units', 'Points', 'Position', [20+depth*40 (counter*70)-20 w-(depth+1)*40 30]);
        hold on;
        set(ahs, 'ytick', []);
        title(names, 'FontUnits', 'Normalized');
        shs = fill([0 0 0 0], [0 0 1 1], [.4 .6 1]);
        ths = text(.05, .5, '0% Complete');
        set(ths, 'FontUnits', 'Normalized');
        ths(1, 2) = text(.45, .5, '');
        set(ths(1, 2), 'FontUnits', 'Normalized');
        axis([0 1 0 1]);
        axis manual;
    else
        depth = depth + 1;
        n = size(names, 1) ;
        for c1 = 0:n-1          
           [timerhs{n-c1, 1} ahs{n-c1, 1} ths{n-c1, 1} shs{n-c1, 1} ps{n-c1, 1} w counter depth] = makebar(names{n-c1, 1}, w, counter, depth);
        end;
        depth = depth - 1;
    end;
end


%--------------------------------------------------------------------
%--------------------------------------------------------------------



function handles = update_report(handles, ps, starter)

% This function generates a progress report.
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'handles':  This input should contain the handles of progress report
            objects as outputted by progress_report().

'ps':       This input is a cell array of the same form as a 'titles' 
            inputted into progress_report(), but with elements of the 
            interval [0, 1] in place of strings, each corresponding 
            to the progress of the matching string. If 'starter' 
            is supplied, this variable instead indicate which 
            progress bar to start or reset. See the code of 
            Assemble() for examples of the form that 'ps' should take 
            in this case.

'starter':  This variable should only be supplied to start or reset a
            progress bar. A value of 1 starts the progress bar then 
            determined by 'ps'; a value of 0 resets that bar. See the 
            code of Assemble() for examples of the form that 'ps' 
            should take in this case.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'handles':  This output is a cell array containing the handles for 
            objects in the progress report.

%--------------------------------------------------------------------
%}

if(exist('starter', 'var'))
    handles{5} = start(handles{5}, ps, starter); 
else

figure(handles{1});
ahs = handles{2};
ths = handles{3};
shs = handles{4};
timerhs = handles{5};
n = size(ps, 1);

for c1 = 0:n-1
    [timerhs{n-c1, 1} ahs{n-c1, 1} ths{n-c1, 1} shs{n-c1, 1} ps{n-c1, 1}] = update(timerhs{n-c1, 1}, ahs{n-c1, 1}, ths{n-c1, 1}, shs{n-c1, 1}, ps{n-c1, 1});    
end;

handles = {handles{1} ahs ths shs timerhs};

end;

end


%--------------------------------------------------------------------
%--------------------------------------------------------------------



function [timerhs, ahs, ths, shs, ps] = update(timerhs, ahs, ths, shs, ps)

% This recursively defined function is for use by the update_report()
% function.

if(isa(ps, 'cell'))
    n = size(ps, 1);
    for c1 = 0:n-1
        [timerhs{n-c1, 1} ahs{n-c1, 1} ths{n-c1, 1} shs{n-c1, 1} ps{n-c1, 1}] = update(timerhs{n-c1, 1}, ahs{n-c1, 1}, ths{n-c1, 1}, shs{n-c1, 1}, ps{n-c1, 1});
    end;
else
    axes(ahs); 
    hold on;
    delete(shs);
    shs = fill([0 ps ps 0], [0 0 1 1], [.4 .6 1]);
    if(isempty(timerhs))
        tstr = '';
    else
        if(ps == 1)
            if(isa(timerhs, 'uint64'))
                timerhs = toc(timerhs);
            end;            
            tstr = ['Ran ' sec2str(timerhs)];           
        else
            if(ps == 0)
                ttime = toc(timerhs);           
                tstr = {[sec2str(ttime) ' Running' ]};
            else
                ttime = toc(timerhs);           
                tstr = {[sec2str(ttime) ' Running' ]; ['~' sec2str(ttime/ps-ttime) ' Remaining']};
            end;
        end;
    end;
    set(ths(1, 1), 'String', [ num2str(floor(1000*ps)/10) '% Complete' ]);
    set(ths(1, 2), 'String', tstr);
    set(ahs, 'Children', [ths(1, 1) ths(1, 2) shs]);
end;



end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function timerhs = start(timerhs, ps, starter)

% This recursively defined function is for use by the update_report()
% function.

if(isa(timerhs, 'cell'))
    timerhs{ps(1, 1), 1} = start(timerhs{ps(1, 1), 1}, ps(1, 2:end), starter);
else
    if(starter)
        timerhs = tic;
    else
        timerhs = [];
    end;
end;
    
    
end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function str = sec2str(sec)

% This function converts a number of seconds into a time string.
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'sec':      This variable should be a double representing a time in
            seconds.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'str':      This output is a string that indicates the amount of time
            respresented by 'sec' in a neater way.

%--------------------------------------------------------------------
%}

if(sec < 60)
    if(abs(ceil(sec)) ~= 1)
        str = [num2str(ceil(sec)) ' Seconds'];
    else
        str = [num2str(ceil(sec)) ' Second'];
    end;
else
    seconds = mod(sec, 60);
    min = (sec-seconds)/60;
    if(min < 60)
        str = [num2str(min) ' Min ' num2str(ceil(seconds)) ' Sec'];
    else
        minutes = mod(min, 60);
        hr = (min - minutes)/60;
        if(hr < 24)
            str = [num2str(hr) ' H ' num2str(minutes) ' M ' num2str(ceil(seconds)) ' S'];
        else
            hours = mod(hr, 24);
            dy = (hr - hours)/24;
            str = [num2str(dy) ' D ' num2str(hours) ' H ' num2str(minutes) ' M ' num2str(ceil(seconds)) ' S'];
        end;
    end;
end;


end


%--------------------------------------------------------------------
%--------------------------------------------------------------------


function errors = check_parameters(parameters, strict)

% This function checks parameters to see if they are appropriate for 
% use as described in [1].
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'parameters':   This variable should be a cell array containing 
            parameters to be checked. See the default values of 
            'parameters' in Solve_Puzzle() (above) for the form this 
            input must take. 

'strict':   This bool should be 1 if empty values in the parameter
            sequence are not to be permitted, and 0 otherwise.

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'pieces':   This output is a n-by-1 cell array each of whose rows is 
            a string containing a message indicating why a particular 
            parameter value is inappropriate for use as described in [1].

%--------------------------------------------------------------------
%}


% Initialize output
errors = {};
if(~(parameters{1} >= 0))
    errors = [errors ; {'alpha must be a nonnegative real number'} ;' '];
end;
if(~(parameters{2} >= 0))
    errors = [errors ; {'beta must be a nonnegative real number'} ;' '];
end;
if(~(parameters{3} >= 0))
    errors = [errors ; {'gamma must be a nonnegative real number'} ;' '];
end;
if(~(parameters{4} > 0))
    errors = [errors ; {'C_1 must be a positive real number'} ;' '];
end;
if(~(parameters{5} > 0))
    errors = [errors ; {'C_2 must be a positive real number'} ;' '];
end;
if(~(parameters{6} > 0))
    errors = [errors ; {'K_1 must be a positive real number'} ;' '];
end;
if(~(parameters{7} > 0))
    errors = [errors ; {'K_2 must be a positive real number'} ;' '];
end;
if(~(parameters{8} > 0))
    errors = [errors ; {'K_3 must be a positive real number'} ;' '];
end;
if(~(parameters{9} > 0))
    errors = [errors ; {'lambda_0 must be a positive real number'} ;' '];
end;
if(~(parameters{10} > 0))
    errors = [errors ; {'lambda_1 must be a positive real number'} ;' '];
end;
if(~(parameters{11} >= 0))
    errors = [errors ; {'nu must be a nonnegative real number'} ;' '];
end;
if(~(parameters{12} > 0))
    errors = [errors ; {'epsilon must be a positive real number'} ;' '];
end;
if(~(parameters{13} >= 0))
    errors = [errors ; {'rho must be a nonnegative real number'} ;' '];
end;
if(~(isa(parameters{14}, 'double') && round(parameters{14}) == parameters{14}  && parameters{14} > 0))
    errors = [errors ; {'j_max must be a positive integer'} ;' '];
end;
if(~(parameters{15} >= 0))
    errors = [errors ; {'eta_1 must be a nonnegative real number'} ;' '];
end;
if(~(parameters{16} >= 0))
    errors = [errors ; {'eta_2 must be a nonnegative real number'} ;' '];
end;
if(~(parameters{17} >= 0))
    errors = [errors ; {'Q_1 must be a nonnegative real number'} ;' '];
end;
if(~(parameters{18} >= 0))
    errors = [errors ; {'Q_2 must be a nonnegative real number'} ;' '];
end;
if(~(parameters{19} >= 0))
    errors = [errors ; {'Q_2^{star} must be a nonnegative real number'} ;' '];
end;
if(~(parameters{20} >= 0))
    errors = [errors ; {'Q_3 must be a nonnegative real number'} ;' '];
end;


parameter_seq = parameters{21};
for c1 = 1:size(parameter_seq, 2)
    if(~(parameter_seq(c1).p_0 >= 0 && parameter_seq(c1).p_0 <= 1) && (~isnan(parameter_seq(c1).p_0) || strict))
        errors = [errors ; {['p_(0,' num2str(c1) ') must belong to the interval [0, 1]']} ;' '];
        parameter_seq(c1).p_0
    end;
    if(~(isa(parameter_seq(c1).m_0, 'double') && round(parameter_seq(c1).m_0) == parameter_seq(c1).m_0  && parameter_seq(c1).m_0 > 0) && (~isnan(parameter_seq(c1).m_0) || strict))
        errors = [errors ; {['m_(0,' num2str(c1) ') must be a positive integer']} ;' '];
    end;
    if(~(parameter_seq(c1).mu_0 >= 0 && parameter_seq(c1).mu_0 <= 1) && (~isnan(parameter_seq(c1).mu_0) || strict))
        errors = [errors ; {['mu_(0,' num2str(c1) ') must belong to the interval [0, 1]']} ;' '];
    end;
    if(~(parameter_seq(c1).K_3 > 0) && (~isnan(parameter_seq(c1).K_3) || strict))
        errors = [errors ; {['K_(3,' num2str(c1) ') must be a positive real number']} ;' '];
    end;
end;


end

%--------------------------------------------------------------------
%--------------------------------------------------------------------


function lengths = vecnorm(points)

% This function calculates the norms of a collection of points in R^2
%{
%--------------------------------------------------------------------
INPUTS
%--------------------------------------------------------------------

'points':   This variable should contain points represented as an 
            n-by-2 matrix, each row of which specifies a point in R^2

%--------------------------------------------------------------------


%--------------------------------------------------------------------
OUTPUTS
%--------------------------------------------------------------------

'lengths':  This output is a n-by-1 matrix whose jth element is the norm of
            the jth point in 'points'

%--------------------------------------------------------------------
%}

lengths = (points(:, 1).^2 + points(:, 2).^2).^(1/2);


end



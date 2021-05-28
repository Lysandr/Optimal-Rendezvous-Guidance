%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new-aiaa.cls, v1.0, 2017/07/31
% Developed by Overleaf
% 
% This class file enables authors to prepare papers
% for submission to AIAA Technical Journals and
% Conferences.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\ProvidesClass{new-aiaa}[2017/07/31, v1.0]

\newif\if@conf
\newif\if@journal
\DeclareOption{conf}{\@conftrue}
\DeclareOption{journal}{\@journaltrue}

\ExecuteOptions{}
\ProcessOptions\relax
\LoadClass[10pt]{article}

% \RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{microtype}
% \RequirePackage{newtxtext,newtxmath}

\RequirePackage[letterpaper,margin=1in]{geometry}
\RequirePackage{enumitem}
\setlist{labelindent=\parindent,leftmargin=*,nosep}
\setlist[enumerate]{label={\arabic*)}}

\RequirePackage[symbol*]{footmisc}

\RequirePackage{setspace}
\if@conf
\fi
\if@journal
  \doublespacing
\fi

\RequirePackage[style]{abstract}
\renewcommand{\abstitlestyle}[1]{}
\renewcommand{\abstracttextfont}{\bfseries\normalsize}
\setlength{\absleftindent}{0.5in}
\setlength{\absrightindent}{0.5in}


\RequirePackage[blocks]{authblk}
\renewcommand{\Authfont}{\normalsize\upshape}
\renewcommand{\Affilfont}{\normalsize\itshape}
\setlength{\affilsep}{0pt}

\renewcommand{\maketitle}{%
  \begingroup
  \renewcommand\and{\\[\baselineskip]}
  \centering\singlespacing
  {\LARGE\bfseries\@title\par}
  \vskip2\baselineskip
  \@author\par
  \endgroup
  \vskip\baselineskip
  \thispagestyle{empty}
  \setcounter{footnote}{0}
}


\RequirePackage[explicit]{titlesec}
\renewcommand{\thesection}{\Roman{section}}
\renewcommand{\thesubsection}{\Alph{subsection}}
\renewcommand{\thesubsubsection}{\arabic{subsubsection}}
\titleformat{\section}
  {\large\bfseries\singlespacing\centering}
  {\thesection.\space}{0pt}{#1}[]
\titlespacing{\section}{0pt}{0.5\baselineskip}{0pt}
\titleformat{\subsection}
  {\normalsize\bfseries\singlespacing}
  {\thesubsection.\space}{0pt}{#1}[]
\titlespacing{\subsection}{0pt}{0.5\baselineskip}{0pt}
\titleformat{\subsubsection}
  {\normalsize\itshape\singlespacing}
  {\thesubsubsection.\space}{0pt}{#1}[]
\titlespacing{\subsubsection}{0pt}{0.5\baselineskip}{0pt}

\RequirePackage{lettrine}
\setlength{\DefaultNindent}{0pt}

\RequirePackage[tableposition=top]{caption}
\renewcommand{\figurename}{Fig.}
\captionsetup*{font=bf,labelsep=quad}
\renewcommand{\arraystretch}{1.2}

\RequirePackage[indentfirst=false,font+=small,leftmargin=0.4in,rightmargin=0pt,vskip=0pt]
{quoting}

\setlength{\abovedisplayskip}{\baselineskip}
\setlength{\belowdisplayskip}{\baselineskip}
\setlength{\abovedisplayshortskip}{\baselineskip}
\setlength{\belowdisplayshortskip}{\baselineskip}

\RequirePackage[sort&compress,numbers]{natbib}
% \bibliographystyle{aiaa}
\bibliographystyle{new-aiaa}
\renewcommand{\bibfont}{\small}

\RequirePackage[hyphens]{url}
\RequirePackage{hyperref}

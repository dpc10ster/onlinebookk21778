\encoding{UTF-8}
\name{RJafroc-package}
\alias{RJafroc-package}
\docType{package}
\title{
  JAFROC analysis for MRMC data
}
\description{
  Software for the analysis of data acquired in observer performance studies conducted using ROC, FROC or ROI multiple reader multiple case (MRMC) data collection paradigms. It is an R implementation of current JAFROC analysis (\url{http://www.devchakraborty.com}) with enhancements including allowing choice between DBM or OR significance testing methods, with Hillis improvements, and choices of several figures of merit, plotting empirical operating characteristics, sample size estimation tools and generating formatted outputs.
}
\details{
  \tabular{ll}{
    Package: \tab RJafroc\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2015-04-09\cr
    License: \tab GPL-3\cr
    URL: \tab \url{http://www.devchakraborty.com}\cr
  }
}

\section{Abbreviations and definitions}{
  \itemize{
      \item treatment/modality: used interchangeably, for example, CT images vs. MRI images of the same patients
      \item reader/observer: used interchangeably, also radiologist
      \item image/case: used interchangeably; a case can consist of several images of the same patient in the same modality
      \item MRMC: multiple reader multiple case (each reader interprets each case in each modality, i.e. fully crossed study design)
      \item DBM: Dorfman-Berbaum-Metz (a significance testing method for detecting a treatment effect in MRMC studies)
      \item DBMH: Hillis modification of DBM
      \item OR: Obuchowski-Rockette (a significance testing method for detecting a treatment effect in MRMC studies)
      \item ORH: Hillis modification of OR
      \item ROC: receiver operating characteristic (a data collection paradigm where each image yields a single rating)
      \item mark: the location of a suspected diseased region 
      \item rating: the level of confidence that the image or the location is diseased, higher numbers indicate increasing confidence in presence of disease
      \item FROC: free-response ROC (a data collection paradigm where each image yields a random number, 0, 1, 2,..., of mark-rating pairs)
      \item AFROC: alternative FROC
      \item	JAFROC: jackknife AFROC: an integrated software suite for analyzing observer performance data
      \item ROI: region-of-interest (each case is divided into a fixed number of region and the reader rates each region)
      \item FOM: figure of merit or quantitative measure of performance
      \item FP: false positive
      \item TP: true positive
      \item FPF: number of FPs divided by number of non-diseased cases
      \item TPF: number of TPs divided by number of diseased cases
      \item SP: specificity, same as \eqn{1-FPF}
      \item SE: sensitivity, same as \eqn{TPF}
      \item ROC operating characteristic: plot of \eqn{TPF} (ordinate) vs. \eqn{FPF}
      \item AUC: trapezoidal area under the ROC curve as estimated by the Wilcoxon statistic
      \item NL: non-lesion localization, of which FP is a special case, i.e., a mark that does not correctly locate any existing localized lesion(s) 
      \item LL: lesion localization, of which TP is a special case, i.e., a mark that correctly locates an existing localized lesion
      \item LLF: number of LLs divided by the total number of lesions
      \item NLF: number of NLs divided by the total number of cases
      \item FROC curve: plot of \eqn{LLF}  (ordinate) vs. \eqn{NLF}
      \item AFROC curve: plot of \eqn{LLF} (ordinate) vs. \eqn{FPF}, where \eqn{FPF} is inferred using highest rating of NL marks on non-diseased cases only
      \item AFROC1 curve: plot of \eqn{LLF} (ordinate) vs. \eqn{FPF1}, where \eqn{FPF1} is inferred using highest rating of NL marks on all cases
      \item JAFROC FOM: trapezoidal area under AFROC curve
      \item JAFROC1 FOM: trapezoidal area under AFROC1 curve
      \item alpha/\eqn{\alpha}: The significance level of the test of the null hypothesis of no treatment effect  
      \item p-value: the probability, under the null hypothesis, that the observed treatment effects, or larger, could occur by chance
      \item NH: The null hypothesis that all treatments effects are zero; rejected if the \eqn{p}-value is smaller than \eqn{\alpha}
      \item RRRC: Analysis that treats both readers and cases as random factors 
      \item RRFC: Analysis that treats readers as random and cases as fixed factors 
      \item FRRC: Analysis that treats readers as fixed and cases as random factors 
      \item ddf: Denominator degrees of freedom of appropriate \eqn{F}-test, the numerator df (ndf) is always number of treatments minus one  
      \item CI: The \eqn{1-\alpha} confidence interval for the stated statistic
      \item I: total number of modalities, indexed by \eqn{i}; \eqn{I} must be at least 2 to perform  null hypothesis testing
      \item J: total number of readers, indexed by \eqn{j}
      \item K1: total number of non-diseased cases, indexed by \eqn{k1}
      \item K2: total number of diseased cases, indexed by \eqn{k2}
      \item K: total number of cases, \eqn{K = K1 + K2}, indexed by \eqn{k}
      \item maxNL:  maximum number of NL marks per case in dataset
      \item maxLL:  maximum number of lesions per case in dataset
    }
}

\section{Dataset}{
  Dataset, an R object, can be created by the user or read from an external data file. \strong{Note}: \bold{the word "dataset" used in this package always 
  represents an R object with following structure}.
  \subsection{Data structure}{
    The dataset is an R list containing 9 elements: \strong{Note}: \code{-Inf} is assigned to any missing/unavailabe element, e.g., an unmarked true lesion.
    \itemize{
      \item \code{NL}: a floating-point array with a dimension of \code{c(I, J, K, maxNL)} that contains the ratings of NL marks for specified modality, reader and case. For ROC datasets FP ratings are assigned to \code{NL} with \code{maxNL = 1}, i.e., the last index is set to 1.
      \item \code{LL}: a floating-point array with a dimension of \code{c(I, J, K2, maxLL)}that contains the ratings of all LL marks for specified modality, reader and case. For ROC datasets TP ratings are assigned to \code{LL} with \code{maxLL = 1}.
      \item \code{lesionNum}: a integer vector with a length of \code{K2}, whose elements indicate the number of lesions in each diseased case. 
      \item \code{lesionID}: a integer array with a dimnsion of \code{c(K2, maxLL)}. \strong{Note} that ratings of lesions in \code{LL} must appear in the same sequence as \code{lesionID} for that case. For example, if the \code{lesionID} field for the first diseased case is \code{c(4, 2, 3, 1)}, i.e., there are 4 lesion on this case labeled 4, 2, 3 and 1, the ratings in \code{LL} for this case must appear in the same sequence, with the first rating corresponding to the lesion labeled 4, the second corresponding to the lesion labeled 2, etc.
      \item \code{lesionWeight}: a floating point array with a dimension of \code{c(K2, maxLL)}, representing the relative importance of detecting each lesion. For each case, the weights must sum to unity. If zero is assigned to all elements of this array, then the software assigns equal weighting, e.g., \code{c(0.5, 0.5)} to an image with two lesions.
      \item \code{maxNL}: the maximum number of \code{NL} marks per case over the entire dataset.
      \item \code{dataType}: a string variable: "ROC", "ROI" or "FROC".
      \item \code{modalityID}: a string vector of length \eqn{I}, which labels the modalities in the dataset.
      \item \code{readerID}: a string vector of length \eqn{J}, which contains the ID of each reader. \strong{Note} that the order of elements in \code{modalityID} and \code{readerID} must match that in \code{NL} and \code{LL}. For example, \code{NL[1, 2, , ]} indicates the ratings of the reader with the second ID in \code{readerID} using the modality with the first ID in \code{modalityID}.
    }
  }
  
  \subsection{Data file format}{
    The package reads JAFROC, MRMC (ROC data only) and iMRMC (ROC data only) data files. The data can be imported by using the function \link{ReadDataFile}. 
    \itemize{
    \item{\strong{JAFROC data file format}}{
    
      The JAFROC data file is an Excel file containing three worksheets (*.xls and *.xlsx are supported): (1) the \samp{Truth} worksheet, (2) the \samp{TP} or lesion localization worksheet and (3) the \samp{FP} or non-lesion localization worksheet. Except for the \samp{Truth} worksheet, where each case must occur at least once, the number of rows in the other worksheets is variable.
      \enumerate{
        \item \samp{Truth} worksheet consists of
          \itemize{
            \item \samp{CaseID}, an integer field uniquely labeling the cases (images). It must occur at least once for each case, and since a case may have multiple lesions, it can occur multiple times, once for each lesion.
            \item \samp{LesionID}, an integer field uniquely labeling the lesions in each case. This field is zero for non-diseased cases.
            \item \samp{Weight}, a floating-point field, which is the relative importance of detecting each lesion. This field is zero for non-diseased cases and for equally weighted lesions; otherwise the weights must sum to unity for each case. Unless a weighted figure of merit is selected, this field is irrelevant.
          }
        \item \samp{TP} worksheet consists of
          \itemize{
              \item \samp{ReaderID}, a string field uniquely labeling the readers (radiologists).
              \item \samp{ModalityID}, a string field uniquely labeling the modalities.
              \item \samp{CaseID}, see \samp{Truth} worksheet. A non-diseased case in this field will generate an error.
              \item \samp{LesionID}, see \samp{Truth} worksheet. An entry in this field that does not appear in the \samp{Truth} worksheet will generate an error. It is the user's responsibility to ensure that the entries in the \samp{Truth} and \samp{TP} worksheets correspond to the same physical lesions.
              \item \samp{TP_Rating}, a positive floating-point field denoting the rating assigned to a particular lesion-localization mark, with higher numbers represent greater confidence that the location is actually a lesion.
            }
        \item \samp{FP} worksheet consists of
          \itemize{
              \item \samp{ReaderID}, see \samp{TP} worksheet.
              \item \samp{ModalityID}, see \samp{TP} worksheet.
              \item \samp{CaseID}, see \samp{TP} worksheet.
              \item \samp{FP_Rating}, a positive floating-point field denoting the rating assigned to a particular non-lesion-localization mark, with higher numbers represent greater confidence that the location is actually a lesion.
            }        
      }
      }
    \item{\strong{MRMC data file format / LABMRMC format}}{
    \itemize{
      \item \emph{Input format for MRMC.} This format is described in the Medical Image Perception Laboratory website, currently \url{http://perception.radiology.uiowa.edu/}.
      \item \emph{LABMRMC data format.} The data file includes following parts. The file must be saved as plain text file with *.lrc extension. All items in the file are separated by one or more blank spaces.
        \enumerate{
          \item The first line is a free text description of the file.
          \item The second line is the name or ID of the first reader.
          \item The third line has the names or IDs of all the modalities. Each name or ID must be enclosed by double quotes(" ").
          \item The fourth line must have the letter (l or s) or word (large or small) for each modality. The letter or word indicates that smaller or larger ratings represent stronger confidence of presence of disease.
          \item The following lines contain the ratings in all modalities, separated by spaces or tabs, of the non-diseased cases, one case per line. The cases must appear in the same order for all readers. Missing value is not allowed. 
          \item After the last non-diseased case insert a line containing the asterisk (*) symbol. 
          \item Repeat steps 5 and 6 for the diseased cases. 
          \item Repeat steps 2, 5, 6 and 7 for the remaining readers. 
          \item The last line of the data file must be a pound symbol (#). 
        }
    
    }
    }
    
    \item{\strong{iMRMC data format}}{
    
      This format is described in the iMRMC website, currently \url{https://code.google.com/p/imrmc/}.
    }
    }
  }
}

\section{Functions}{
  \itemize{    
    \item \code{\link{DBMHAnalysis}}: Performs Dorfman-Berbaum-Metz analysis with Hillis improvements for the specified dataset.
    \item \code{\link{EmpiricalOpCharac}}: Plot empirical curves for specified modalities and readers in the dataset.
    \item \code{\link{FigureOfMerit}}: Calculate the figure of merit for each reader using each modality.
    \item \code{\link{FROC2HrROC}}: Convert an FROC dataset to a highest rating inferred ROC dataset.
    \item \code{\link{ORHAnalysis}}: Performs Obuchowski-Rockette analysis with Hillis improvements for the specified dataset.
    \item \code{\link{OutputReport}}: Save the results of the analysis to a text file.  
    \item \code{\link{PowerGivenJK}}: Calculate the statistical power with the given number of readers, number of cases and DBM or OR variances components.
    \item \code{\link{PowerTable}}: Calculate required sample size for the specified dataset with given significance level, effect size and desired power.
    \item \code{\link{ReadDataFile}}: Read the dataset that will be analysis from data file.
    \item \code{\link{SampleSizeGivenJ}}: Calculate required number of cases with the given number of readers and DBM variances components.
    \item \code{\link{SaveDataFile}}: Save data file in specified format.
  }
}

\author{
  Xuetong Zhai, Dev Chakraborty

  Maintainer: Xuetong Zhai <xuetong.zhai@gmail.com>
}
\references{
\bold{Basics of ROC}

  Metz, C. E. (1978). Basic principles of ROC analysis. In Seminars in nuclear medicine (Vol. 8, pp. 283–298). Elsevier. 
  
  Metz, C. E. (1986). ROC Methodology in Radiologic Imaging. Investigative Radiology, 21(9), 720. 
  
  Metz, C. E. (1989). Some practical issues of experimental design and data analysis in radiological ROC studies. Investigative Radiology, 24(3), 234.
  
  Metz, C. E. (2008). ROC analysis in medical imaging: a tutorial review of the literature. Radiological Physics and Technology, 1(1), 2–12.
  
  Wagner, R. F., Beiden, S. V, Campbell, G., Metz, C. E., & Sacks, W. M. (2002). Assessment of medical imaging and computer-assist systems: lessons from recent experience. Academic Radiology, 9(11), 1264–77. 
  
  Wagner, R. F., Metz, C. E., & Campbell, G. (2007). Assessment of medical imaging systems and computer aids: a tutorial review. Academic Radiology, 14(6), 723–48. 

\bold{DBM/OR methods and extensions}

  DORFMAN, D. D., BERBAUM, K. S., & Metz, C. E. (1992). Receiver operating characteristic rating analysis: generalization to the population of readers and patients with the jackknife method. Investigative Radiology, 27(9), 723.
  
  Obuchowski, N. A., & Rockette, H. E. (1994). HYPOTHESIS TESTING OF DIAGNOSTIC ACCURACY FOR MULTIPLE READERS AND MULTIPLE TESTS: AN ANOVA APPROACH WITH DEPENDENT OBSERVATIONS. Communications in Statistics-Simulation and Computation, 24(2), 285–308. 

  Hillis, S. L., Berbaum, K. S., & Metz, C. E. (2008). Recent developments in the Dorfman-Berbaum-Metz procedure for multireader ROC study analysis. Academic Radiology, 15(5), 647–61.
  
  Hillis, S. L., Obuchowski, N. A., & Berbaum, K. S. (2011). Power Estimation for Multireader ROC Methods: An Updated and Unified Approach. Acad Radiol, 18, 129–142.
  
  Hillis, S. L. S. L. (2007). A comparison of denominator degrees of freedom methods for multiple observer ROC analysis. Statistics in Medicine, 26(3), 596–619. 

\bold{FROC paradigm}

  Chakraborty, D. P., & Berbaum, K. S. (2004). Observer studies involving detection and localization: modeling, analysis, and validation. Medical Physics, 31(8), 1–18. 
  
  Chakraborty, D. P. (2006). A search model and figure of merit for observer data acquired according to the free-response paradigm. Physics in Medicine and Biology, 51(14), 3449–62. 
  
  Chakraborty, D. P. (2006). ROC curves predicted by a model of visual search. Physics in Medicine and Biology, 51(14), 3463–82.
  
  Chakraborty, D. P. (2011). New Developments in Observer Performance Methodology in Medical Imaging. Seminars in Nuclear Medicine, 41(6), 401–418.
  
  Chakraborty, D. P. (2013). A Brief History of Free-Response Receiver Operating Characteristic Paradigm Data Analysis. Academic Radiology, 20(7), 915–919. 
  
  Chakraborty, D. P., & Yoon, H.-J. (2008). Operating characteristics predicted by models for diagnostic tasks involving lesion localization. Medical Physics, 35(2), 435. 

\bold{ROI paradigm}

  Obuchowski, N. A., Lieber, M. L., & Powell, K. A. (2000). Data analysis for detection and localization of multiple abnormalities with application to mammography. Academic Radiology, 7(7), 553–4; discussion 554–6. 
}

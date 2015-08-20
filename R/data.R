#' Krackhardt Office Cognitive Social Structure
#'
#' David Krackhardt collected cognitive social structure data
#' from 21 management personnel in a high-tech, machine manufacturing
#' firm to assess the effects of a recent management intervention program.
#' The relation queried was “Who is a friend of X?" (krackfr). Each person indicated
#' not only his or her own advice and friendship relationships, but also
#' the relations he or she perceived among all other managers, generating a
#' full 21 by 21 matrix of adjacency ratings from each person in the group.
#'
#' @format An array with 21 matrices which are each 21 by 21:
#' \describe{
#'   \item{kfr}{array}
#' }
#' @source Krackhardt D. (1987). Cognitive social structures. Social Networks, 9, 104-134.
#' @examples
#' data(kfr)
#' dim(kfr)
"kfr"

#' Sampson–Monastery
#'
#' Sampson recorded the social interactions among a group of monks
#' while resident as an experimenter on vision, and collected numerous
#' sociometric rankings. During his stay, a political “crisis in the cloister"
#' resulted in the expulsion of four monks (Nos. 2, 3, 17, and 18) and the
#' voluntary departure of several others - most immediately,
#' Nos. 1, 7, 14, 15, and 16. (In the end, only 5, 6, 9, and 11 remained).
#'
#' Most of the present data are retrospective, collected after the breakup
#' occurred. They concern a period during which a new cohort entered the
#' monastery near the end of the study but before the major conflict began. T
#' he exceptions are "liking" data gathered at three times: LikeT1 to LikeT3 -
#' that reflect changes in group sentiment over time (LikeT3 was
#' collected in the same wave as the data described below).
#' Information about the senior monks was not included.
#'
#' Four relations are coded, with separate matrices for
#' positive and negative ties on the relation. Each member ranked
#' only his top three choices on that tie.
#' The relations are esteem (Esteem) and disesteem (Disesteem),
#' liking (LikeT*) and disliking (Dislike), positive influence (Influence)
#'  and negative influence (NegInfluence), praise (Praise) and blame (Blame).
#'  In all rankings 3 indicates the highest or first choice and 1 the last
#'  choice. (Some subjects offered tied ranks for their top four choices).
#'
#' @format An list of 18 by 18 valued matrices:
#' \describe{
#' \item{ Esteem }{matrix}
#' \item{ Disesteem }{matrix}
#' \item{ Influence }{matrix}
#' \item{ NegInfluence }{matrix}
#' \item{ LikeT1 }{matrix}
#' \item{ LikeT2 }{matrix}
#' \item{ LikeT3 }{matrix}
#' \item{ Dislike }{matrix}
#' \item{ Praise }{matrix}
#' \item{ Blame }{matrix}
#' }
#' @source Sampson, S. (1969). Crisis in a cloister. Unpublished doctoral dissertation, Cornell University.
#' @references Breiger R., Boorman S. and Arabie P. (1975). An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. Journal of Mathematical Psychology, 12, 328-383.
#' @examples
#' data(sampson)
#' length(sampson)
"sampson"

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

#' Coleman's High School Friendship Data
#'
#' ames Coleman (1964) reports research on self-reported friendship ties
#' among 73 boys in a small high school in Illinois over the 1957-1958
#' academic year. Networks of reported ties for all 73 informants are
#' provided for two time points (fall and spring).
#'
#' Both networks reflect answers to the question,
#' “What fellows here in school do you go around with most often?”
#' with the presence of an (i,j,k) edge indicating that j nominated
#' k in time period i.
#'  The data are unvalued and directed; although the self-reported ties
#'  are highly reciprocal, unreciprocated nominations are possible.
#'  It should be noted that, although this data is usually described as
#'  “friendship,” the sociometric item employed might be more accurately
#'  characterized as eliciting “frequent elective interaction.” This should
#'  be borne in mind when interpreting this data.
#'
#' @format An array with 2 matrices which are each 73 by 73:
#' \describe{
#'   \item{coleman}{array}
#' }
#' @source Coleman, J. S. (1964). Introduction to Mathermatical Sociology. New York: Free Press.
#' @examples
#' data(coleman)
#' dim(coleman)
#' dimnames(coleman)[1]
"coleman"

#' Kapferer–Tailor Shop
#'
#' Bruce Kapferer (1972) observed interactions in a tailor shop in
#' Zambia (then Northern Rhodesia) over a period of ten months.
#' His focus was the changing patterns of alliance among workers
#' during extended negotiations for higher wages.
#'
#' @format A network object:
#' \describe{
#'   \item{kaptail.ins}{network}
#' }
#' @source Kapferer B. (1972). Strategy and transaction in an African factory. Manchester: Manchester University Press.
#' @examples
#' data(kaptail.ins)
#' kaptail.ins
"kaptail.ins"

#' Thurman Office Network
#'
#' Thurman spent 16 months observing the interactions among
#' employees in the overseas office of a large international
#' corporation. During this time, two major disputes erupted
#' in a subgroup of fifteen people. Thurman analyzed the
#' outcome of these disputes in terms of the network of
#' formal and informal associations among those involved.
#'
#' @format A network object:
#' \describe{
#'   \item{thuroff.int}{network}
#' }
#' @source Thurman B. (1979). In the office: Networks and coalitions. Social Networks, 2, 47-63.
#' data(thuroff.int)
#' thuroff.int
"thuroff.int"

#' Correlates of War 1993: Militarized interstate disputes
#'
#'
#'
#' @format A network object:
#' \describe{
#'   \item{mids_1993}{network}
#' }
#' @source  \url{http://www.correlatesofwar.org/}.
#' data(mids_1993)
#' mids_1993
"mids_1993"

#' Correlates of War 1993: Contiguity among nations
#'
#'
#'
#' @format A network object:
#' \describe{
#'   \item{contig_1993}{network}
#' }
#' @source  \url{http://www.correlatesofwar.org/}.
#' data(contig_1993)
#' contig_1993
"contig_1993"

#' Correlates of War 1993: Alliances among nations
#'
#'
#'
#' @format A network object:
#' \describe{
#'   \item{alliances_1993}{network}
#' }
#' @source  \url{http://www.correlatesofwar.org/}.
#' data(alliances_1993)
#' alliances_1993
"alliances_1993"

#' Correlates of War 1993: Contiguity among nations
#'
#'
#'
#' @format A network object:
#' \describe{
#'   \item{contig_1993}{network}
#' }
#' @source  \url{http://www.correlatesofwar.org/}.
#' data(contig_1993)
#' contig_1993
"contig_1993"

#' Southern Women
#'
#' This is the incidence matrix for the famous ``Southern Women" data set from Davis,
#' Gardner, and Gardner's Deep South study of class and social interaction in
#' Natchez, MS during the mid-1930s. The matrix shows the attendance of 18 women at
#' 14 informal social events during a nine-month observation period, based on
#' ``interviews, the records of participant observers, guest lists, and the newspapers" (DGG, p. 149).
#'
#' @format A matrix object:
#' \describe{
#'   \item{sw.incidence}{matrix}
#' }
#' @source Davis, A., Gardner, B. B., and Gardner, M. R. (2009). Deep South: A social anthropological study of caste and class. Univ of South Carolina Press.
#' data(sw.incidence)
#' sw.incidence
"sw.incidence"


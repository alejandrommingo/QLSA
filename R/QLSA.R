#' gallitoContour
#'
#' gallitoContour is a function that will extract a word contour from Gallito API LSA semantic space. It takes
#' a word or a tupple (separated by "_") and returns a data frame with the contour of that word.
#' @param word A string or a tupple of strings separated by "_" indicating the word for which you want to
#' extract the contour
#' @param neighbors The number of neighbors inside the contour of the word. By default `neighbors = 100`.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space
#' @return  A data frame with the contour of the word is returned
#' @export
gallitoContour = function(word, gallitoCode, neighbors = 100){

  k <- 300 # Define K dimensions
  n = neighbors # Define N neighbors
  matriz <- matrix(0, n, k)


  # Take each word in the file
  text = word

  # Prepare the Gallito API query
  text2 <- "<getNearestNeighboursList xmlns='http://tempuri.org/'><code>"
  text2 <- paste(text2, gallitoCode, sep = "")
  text2 <- paste(text2, "</code><a>", sep = "")
  text2 <- paste(text2, text, sep = "")
  text2 <- paste(text2, "</a><txtNegate></txtNegate><n>", sep = " ")
  text2 <- paste(text2, as.character(n), sep = "")
  text2 <- paste(text2, "</n><lenght_biased>false</lenght_biased></getNearestNeighboursList>", sep = "")
  text2 <- enc2utf8(text2)

  # Make te query
  cadena <- httr::POST("http://psicoee.uned.es/quantumlikespace_spanish/Service.svc/webHttp/getNearestNeighboursList", body = text2, httr::content_type("text/xml"))

  # Correct the chain of characters
  txt <- gsub("&lt;", "<", cadena)
  txt <- gsub("&gt;", ">", txt)
  txt <- gsub(",", ".", txt)

  # Parse the text to xml
  doc3 = XML::xmlTreeParse(txt, useInternal = TRUE, encoding = "UTF-8")

  # Take each vector of neighbors from Gallito API
  vector <- XML::getNodeSet(doc3, "//r:term/text()", c(r = "http://tempuri.org/"))
  text <- XML::xmlValue(vector[[1]])
  text <- enc2utf8(text)
  words = c()

  for (j in 1:length(vector)){

    # Find each neighbor vector in the semantic space

    # Take each neighbor
    text <- XML::xmlValue(vector[[j]])
    text <- enc2utf8(text)

    # Prepare the Gallito API Query
    text2 <- "http://psicoee.uned.es/quantumlikespace_spanish/Service.svc/webHttp/getVectorOfTerm?code="
    text2 <- paste(text2, gallitoCode, sep = "")
    text2 <- paste(text2, "&a=", sep = "")
    text2 <- paste(text2, text, sep = "")
    text2 <- enc2utf8(text2)

    # Make the query
    cadena <- httr::GET(text2, httr::content_type("text/xml"))

    # Correct the chain of characters
    txt <- gsub("&lt;", "<", cadena)
    txt <- gsub("&gt;", ">", txt)
    txt <- gsub(",", ".", txt)

    # Parse the text to xml
    doc3 = XML::xmlTreeParse(txt, useInternal = TRUE)

    # Extract the vector from the LSA
    valueByText = as.numeric(XML::xpathApply(doc3, "//r:dim/text()", XML::xmlValue, namespaces = c(r = "http://schemas.microsoft.com/2003/10/Serialization/")))

    # Store all the vectors in a matrix to create the contour
    matriz[j,] = valueByText
    words = cbind(words, XML::xmlValue(vector[[j]]))

  }

  # Using a data frame instead of a matrix
  df = as.data.frame(matriz)
  rownames(df) = words
  return(df)

}


#' wordVector
#'
#' wordVector is a function that will extract a word vector from Gallito API LSA semantic space. It takes a word or a tupple (separated by "_") and returns the vector of that word in the LSA semanti space.
#' @param word A string or a tupple of strings separated by "_" indicating the word you want to extract.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space
#' @return The word vector is returned
#' @export
wordVector = function(word, gallitoCode){
  query <- "http://psicoee.uned.es/quantumlikespace_spanish/Service.svc/webHttp/getVectorOfTerm?code="
  query <- paste(query, gallitoCode, sep = "")
  query <- paste(query, "&a=", sep = "")
  query <- paste(query, word, sep = "")
  query <- enc2utf8(query)
  cadena <- httr::GET(query, httr::content_type("text/xml"))
  txt <- gsub("&lt;", "<", cadena)
  txt <- gsub("&gt;", ">", txt)
  txt <- gsub(",", ".", txt)
  doc3 = XML::xmlTreeParse(txt, useInternal = TRUE)
  valueByText = as.numeric(XML::xpathApply(doc3, "//r:dim/text()", XML::xmlValue, namespaces = c(r = "http://schemas.microsoft.com/2003/10/Serialization/")))
  return(valueByText)
}

#' subspaceGeneration
#'
#' subspaceGeneration is a function that decides the number of dimensions a contour deserves using the `paran` package,
#' compute an EFA (Exploratory Factor Analysis) with that number of dimensions using the `psych` package,
#' and reorthogonalize the solution to define a basis for the new subspace using the `pracma` package.
#' @param word A string or a tupple of strings separated by "_" indicating the word for which you want to
#' define the subspace.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space.
#' @param min_cosine The minimum cosine the function will use to return the plots of similar words for
#' each new dimension of the subspace. By default `min_cosine = 0.5`.
#' @param min_reilability The minimum reilability the function will consider to decide
#' that a reorthogonalized dimension is equivalent to the dimension extracted in the factorial solution.
#' By default `min_reilability = 0.85`.
#' @return The function will return a list with the subspace as `subspace`, the reilability test as `reilability_test`,
#' the subspace graphical information as `subspace_info` and the EFA results as `EFA_info`.
#' @export
subspaceGeneration = function(word, gallitoCode, min_cosine = 0.5, min_reilability = 0.85){

  # Extract the contour from Gallito API
  word_contour = QLSA::gallitoContour(word, gallitoCode)

  # Perform parallel analysis and store the deserved dimensionality
  dim = paran::paran(t(word_contour), quietly = TRUE, status = FALSE)$Retained

  # Perform EFA with that number of dimensions
  word_afe = psych::fa(t(word_contour), nfactors = dim, rotate = "oblimin", fm = "ml")

  # Reorthogonalize the solution using the "subspaceBasis" function
  word_subspace = QLSA::subspaceBasis(fact_anal = word_afe, word_contour = word_contour, n_dim = dim, n_GS_rotations = 100)
  word_subspace_reilability = word_subspace$reilability_test

  # Test the reilability: While there is a value less than 0.85 we keep reducing the number of dimensions until all dimensions are deserved
  while (min(word_subspace_reilability) < min_reilability) {
    dim = dim - 1
    word_afe = psych::fa(t(word_contour), nfactors = dim, rotate = "oblimin", fm = "ml")
    word_subspace = QLSA::subspaceBasis(fact_anal = word_afe, word_contour = word_contour, n_dim = dim, n_GS_rotations = 100)
    word_subspace_reilability = word_subspace$reilability_test
  }

  # Isolate the Subspace
  word_SP = t(word_subspace$orthogonal_subspace)

  # Extract Subspace info using the "subspaceInfo" function (wordclouds, barplots and dimensions terms)
  SP_info = QLSA::subspaceInfo(word_contour, word_SP, min_cosine = min_cosine)

  # Final list to return
  subspaceList = list(subspace = word_SP, reilability_test = word_subspace_reilability, subspace_info = SP_info, EFA_info = word_afe)
  return(subspaceList)
}


#' unidimensionalProjector
#'
#' unidimensionalProjector is a function that takes a unidimensional subspace and calculates its projector using
#' the outter product operation. You can use the `wordVector` function to extract a word vector and calculate its
#' unidimensional projector.
#' @param vector A unidimensional subspace in a vector format.
#' @return The function will return a matrix with the projector of the subspace.
#' @export
unidimensionalProjector = function(vector){
  # Vector outer product
  vector %*% t(vector)
}

#' multidimensionalProjector
#'
#' multidimensionalProjector is a function that takes a multidimensional subspace and calculates its projector using
#' the sum of the outter products of its vectors. You can use the `subspaceGeneration` function from `QLSA` library to
#' define a subspace in order to calculate its multidimensional projector
#' @param subspace A multidimensional subspace in a matrix format.
#' @return The function will return a matrix with the projector of the subspace.
#' @export
multidimensionalProjector = function(subspace){
  x = as.matrix(subspace)
  P = matrix(0,nrow = ncol(x),ncol = ncol(x))
  # For each dimension, calculate outter product and sum all
  for (i in 1:nrow(x)){
    P = P + (x[i,] %*% t(x[i,]))
  }
  return(P)
}

#' contextualDistance
#'
#' contextualDistance is a function that takes two words projectors and a state vector and returns the contextual
#' distance between them. To see more details about the contextual distance go to Gabora and Aerts (2002).
#' @param word_a The first word the function will evaluate.
#' @param word_b The second word the function will evaluate.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space.
#' @param neutral_state By default the function will take the neutral state between the two words subspaces
#' to evaluate the similarity
#' @param state The initial state the function will use to calculate the similarity between the two words in case
#' `neutral_state = FALSE`.
#' @return The function will return a value between 0 and +Inf indicating the distance between the two words.
#' @export
contextualDistance = function(word_a, word_b, gallitoCode, neutral_state = TRUE, state){

  word_a_sbs = QLSA::subspaceGeneration(word_a, gallitoCode)$subspace
  word_a_PR = QLSA::multidimensionalProjector(word_a_sbs)
  word_b_sbs = QLSA::subspaceGeneration(word_b, gallitoCode)$subspace
  word_b_PR = QLSA::multidimensionalProjector(word_b_sbs)

  if (neutral_state == TRUE){
    state = QLSA::neutralState(word_a, word_b, gallitoCode)
    return(sqrt(2*(1-sqrt((norm(word_b_PR%*%word_a_PR%*%state, type = "2"))^2))))
  }
  else{
    return(sqrt(2*(1-sqrt((norm(word_b_PR%*%word_a_PR%*%state, type = "2"))^2))))
  }


}

#' quantumSimilarity
#'
#' quantumSimilarity is a function that takes two words projectors and a state vector and returns the quantum
#' similarity between them. To see more details about the quantum similarity go to Pothos and Busemeyer (2011).
#' @param word_a The first word the function will evaluate.
#' @param word_b The second word the function will evaluate.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space.
#' @param neutral_state By default the function will take the neutral state between the two words subspaces
#' to evaluate the similarity
#' @param state The initial state the function will use to calculate the similarity between the two words in case
#' `neutral_state = FALSE`.
#' @return The function will return a value between 0 and 1 indicating the similarity between the two words.
#' @export
quantumSimilarity = function(word_a, word_b, gallitoCode, neutral_state = TRUE, state){

  word_a_sbs = QLSA::subspaceGeneration(word_a, gallitoCode)$subspace
  word_a_PR = QLSA::multidimensionalProjector(word_a_sbs)
  word_b_sbs = QLSA::subspaceGeneration(word_b, gallitoCode)$subspace
  word_b_PR = QLSA::multidimensionalProjector(word_b_sbs)

  if (neutral_state == TRUE){
    state = QLSA::neutralState(word_a, word_b, gallitoCode)
    return(norm(word_b_PR %*% word_a_PR %*% state ,type = "2")^2)
  }
  else{
    return(norm(word_b_PR %*% word_a_PR %*% state ,type = "2")^2)
  }

}

#' neutralState
#'
#' neutralState is a function that takes two words projectors and estimates the neutral state
#' vector in between. To see more details about the quantum similarity go to Pothos and Busemeyer (2011).
#' @param word_a The first word the function will evaluate.
#' @param word_b The second word the function will evaluate.
#' @param gallitoCode Gallito API password to extract information from the LSA semantic space.
#' @return The function will return a value between 0 and 1 indicating the similarity between the two words.
#' @export
neutralState = function(word_a, word_b, gallitoCode){

  word_a_sbs = QLSA::subspaceGeneration(word_a, gallitoCode)$subspace
  word_a_PR = QLSA::multidimensionalProjector(word_a_sbs)
  word_b_sbs = QLSA::subspaceGeneration(word_b, gallitoCode)$subspace
  word_b_PR = QLSA::multidimensionalProjector(word_b_sbs)

  word_a_vector = QLSA::wordVector(word_a, gallitoCode)
  word_b_vector = QLSA::wordVector(word_b, gallitoCode)

  # Estimate Neutral State
  intermediate_vector = (word_a_vector+word_b_vector)/(sqrt(sum((word_a_vector+word_b_vector)^2)))
  stVecOptim = optim(par = intermediate_vector, fn = stateOptim, PA=word_a_PR, PB=word_b_PR)
  state_vector = stVecOptim$par/sqrt(sum(stVecOptim$par^2))

  return(state_vector)
}







# Background functions.

subspaceBasis = function(fact_anal, word_contour, n_dim, n_GS_rotations){

  # Extract the factor scores
  subspace_fa = psych::factor.scores(f = fact_anal, x = t(word_contour))

  # Create a orthogonalizations list to store all the reorthogonalized matrices
  subspace = subspace_fa$scores
  colnames(subspace) = 1:n_dim
  orthogonalizations_list = list()

  # Re-orthogonalize the subspace using Gram Schmidt several times
  for (j in 1:n_GS_rotations){
    temp_matrix = subspace[,sample(ncol(subspace))]
    temp_orthogonalized_matrix = pracma::gramSchmidt(temp_matrix)$Q
    colnames(temp_orthogonalized_matrix) = colnames(temp_matrix)
    orthogonalized_matrix = temp_orthogonalized_matrix[,order(colnames(temp_orthogonalized_matrix))]
    orthogonalizations_list[[j]] = orthogonalized_matrix
  }

  final_orthogonalized_matrix = matrix(0L, nrow = 300, ncol = n_dim)

  # Obtain the average reorthogonalized matrix
  for (j in 1:length(orthogonalizations_list)){
    final_orthogonalized_matrix = final_orthogonalized_matrix + orthogonalizations_list[[j]]
  }
  final_orthogonalized_matrix = final_orthogonalized_matrix/n_GS_rotations
  final_orthogonalized_matrix = pracma::gramSchmidt(final_orthogonalized_matrix)$Q

  # Obtain the correlations between the original subspace and the reorthogonalized one to test the reilability
  GS_reilability_test = diag(cor(final_orthogonalized_matrix,subspace))

  # Return the final reorthogonalized matrix and the reilability test
  return(list(orthogonal_subspace = final_orthogonalized_matrix, reilability_test = GS_reilability_test))
}

subspaceInfo = function(word_contour, subspace, min_cosine = 0.5){

  # Create a list to store terms for each dimensions
  dim_terms = list()

  # Iterate through all dimensions of a subspace to extract the most similar terms
  for (i in 1:nrow(subspace)){
    count = 0
    df = data.frame(matrix(0, ncol = 3))
    colnames(df) = c("word","cosine","angle")

    for (j in 1:nrow(word_contour)){
      cos = lsa::cosine(as.vector(t(word_contour[j,])), subspace[i,])
      if (cos > min_cosine){
        count = count + 1
        df[count,1] = rownames(word_contour)[j]
        df[count,2] = cos
      }
    }

    # Store all terms for each dimension in a data frame
    dim_terms[[i]] = df[order(df$cosine, decreasing = T),]
  }

  # Create a list to store the wordclouds
  wordclouds = list()

  # Iterate through dimensions to create wordclouds
  for (i in 1:nrow(subspace)){
    df = dim_terms[[i]]
    df[,3] = 90 * sample(c(0, 1), nrow(df), replace = TRUE, prob = c(60, 40))
    word_plot = ggplot2::ggplot(df, ggplot2::aes(label = word, size = cosine, angle = angle)) +
      ggwordcloud::geom_text_wordcloud_area() +
      ggplot2::theme_minimal()
    wordclouds[[i]] = word_plot
  }

  # Create a list for barplots
  barplots = list()

  # Iterate through dimensions to create barplots
  for (i in 1:nrow(subspace)){
    df = dim_terms[[i]]
    df$word <- factor(df$word, levels = df$word[order(df$cosine)])
    bar_plot = ggplot2::ggplot(df, ggplot2::aes(x = word, y = cosine)) +
      ggplot2::geom_bar(stat="identity")+
      ggplot2::coord_flip()
    barplots[[i]] = bar_plot
  }

  # Return dimension terms, wordclouds and barplots
  return(list(terms_list = dim_terms, wordclouds = wordclouds, barplots = barplots))
}

stateOptim = function(state,PA,PB){
  abs(norm(PA%*% state, "2")^2 - norm(PB%*% state, "2")^2)
}

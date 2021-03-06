% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XGSA.R
\name{generate_sparse_matrix}
\alias{generate_sparse_matrix}
\title{generate_sparse_matrix}
\usage{
generate_sparse_matrix(homology_table)
}
\arguments{
\item{homology_table}{Table of homology mapping between two species, as retrieved by the function get_homology_matrix()}
}
\value{
This function return a sparse matrix
}
\description{
This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
It uses the simple_triplet_matrix structure from the package 'slam'.
}
\details{
This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
It uses the simple_triplet_matrix structure from the package 'slam'.
}
\examples{
Not used trivially

}


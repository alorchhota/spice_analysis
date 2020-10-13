def get_matrix_size_dependent_default_thread(wildcards):
  ng = int(wildcards.n_genes)
  if ng <= 1500:
    th = 1
  elif ng <= 5000:
    th = 2
  else:
    th = 1
  return th

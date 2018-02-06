subroutine wkmeans(X, k, W, iter_max, n, m, centroids, clusters, ss)
  ! input (and output variables)
  integer, intent(in) :: n, m, k, iter_max
  double precision, intent(inout) :: X(n,m), W(n)
  double precision, intent(inout) :: centroids(k,m)
  integer, intent(inout) :: clusters(n)
  double precision, intent(inout) :: ss

  ! local work variables
  integer :: i_n, i_m, i_k, iter
  logical :: idx_k(n)
  double precision :: this_sqdist(k), sqdist(n), new_centroids(k,m)

! ALGORITHM
! start from k provided centroids
! iterate
!   compute the squared Euclidean distance between each of the k centroids and each of the n points in X
!   for each of the n points in X, determine the nearest centroid
!   compute new centroids as the weighted barycenter of the points
! until the centroids do not move
! compute the sum of squared residual distances (between each centroid and the points in this cluster)
!
! Optionally, repeat that r times and keep the solution with the smalled sum of squares

  do iter = 1, iter_max, 1
    ! Assign a cluster number to each input point
    ! for each point
    do i_n = 1, n, 1
      ! compute  distance between this point and each centroid
      do i_k = 1, k, 1
        this_sqdist(i_k) = sum((X(i_n,:) - centroids(i_k,:))**2)
      end do
      ! find the nearest centroid
      clusters(i_n) = minloc(this_sqdist, 1)
      ! and store the distance from the current point to this centroid
      sqdist(i_n) = this_sqdist(clusters(i_n))
    end do

    ! Compute new centroids
    ! for each cluster
    do i_k = 1, k, 1
      ! find which points belong to this cluster
      idx_k = (clusters == i_k)
      ! compute the weighted barycenter in each dimension
      do i_m = 1, m, 1
        new_centroids(i_k, i_m) = sum(X(:, i_m)*W, mask=idx_k) / sum(W, mask=idx_k)
      end do
    end do
    
    ! Stop when centroids converge
    if ( abs(sum(new_centroids - centroids)) < 10.**(-6) ) then
      ! TODO this absolute limit is not a good idea since it should depend on the actual scale of the data
      centroids = new_centroids
      exit
    else
      centroids = new_centroids
    end if
    
  end do

  ! store weighted total squared distances to centroids
  ss = sum(sqdist*W)
    
end subroutine wkmeans

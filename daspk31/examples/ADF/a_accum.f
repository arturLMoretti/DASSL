      subroutine a_accum_d_0(l,lhs)
      integer l
      double precision lhs(l)
      do i = 1, l
      lhs(i) = 0.0
      enddo
      end

      subroutine a_accum_d_2(l,lhs,w_1,v_1)
      integer l
      double precision lhs(l)
      double precision v_1(l)
      double precision w_1
      do i = 1, l
      lhs(i) = w_1*v_1(i)
      enddo
      end

      subroutine a_accum_d_4(l,lhs,w_1,v_1,w_2,v_2)
      integer l
      double precision lhs(l)
      double precision v_1(l)
      double precision w_1
      double precision v_2(l)
      double precision w_2
      do i = 1, l
      lhs(i) = w_1*v_1(i) + w_2*v_2(i)
      enddo
      end


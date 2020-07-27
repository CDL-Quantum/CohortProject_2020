      subroutine dmatvecprod(prod,arg1,arg2,size)
      implicit none
      integer size,i
      real*8 prod(size),arg1(size),arg2(size)

      do 10 i=1,size
         prod(i)=arg1(i)*arg2(i)
 10   continue
      return
      end
      subroutine dotprod(arg1,arg2,size,dot)
      implicit none
      integer size,i
      real*8 dot,arg1(size),arg2(size),ddot
c      dot=ddot(size,arg1,1,arg2,1)
      dot=0.d0
      do 10 i=1,size
         dot=dot+arg1(i)*arg2(i)
 10   continue
      return
      end
      subroutine cdotprod(arg1,arg2,size,dot)
      implicit none
      integer size,i
      complex*16 dot,arg1(size),arg2(size),zdotu
c      dot=zdotu(size,arg1,1,arg2,1)
      dot=(0.d0,0.d0)
      do 10 i=1,size
         dot=dot+dconjg(arg1(i))*arg2(i)
 10   continue
      return
      end

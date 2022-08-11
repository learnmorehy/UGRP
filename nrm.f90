	program nrm
	
	implicit none
	
	real(8) :: pi, s, m, n, alpha, alpha1, angle, m_t
	real(8) :: r, x, y, dx, dy 
	real(8) :: x1, x2, y1, y2, m1, m2, m3
	real(8) :: theta1, theta2, f1, f2, e1, e2, a, b, c, t1, t2, k
	real(8) :: s_psos, p_psos, arc, arc1, arc2, arc3, h, theta_arc
	real(8) :: params(11), g1, g2, g3, e, s_num, p_num, num, n_i, n_t
	real(8) :: k1, k2, x_i, y_i, v1(2), v2(2), v3(2), v4(2)
	integer :: i0, i1, i2, i3, i4, i5, bs, i, j, row, col, tmp_r, tmp_c
	real(8) :: r_s, r_p, ref, emi, theta_t, theta_i, intensity, tmp_s, tmp_p, tmp_x, tmp_y
	real(8) :: l, angle_w, angle_f, mat2(370), m_f, v_f(2), q, a_f, b_f, c_f, e_f1, e_f2, x_f, y_f, test
	real(8), DIMENSION(:,:), allocatable :: mat
	character(20) :: dummy(11)
	common /coeff/ g1,g2,g3, l
	
	
	open(unit = 9873, file = "./params.dat", status = "OLD")
	open(unit = 9,file = "./ray.dat")
	open(unit = 1111,file = "./psos.dat")
	open(unit = 1234, file = "./s.dat")
	open(unit = 4321, file = "./spd.dat")
	open(unit = 123, file = "./ffield.dat")
	open(unit = 321, file = "./test_f.dat")
	11 format(X,4E)
    
    pi = 2d0*dasin(1d0)       
	
	do i0 = 1, 11       
		read(9873,*) dummy(i0), params(i0)       
	end do
	
	g1 = params(1)
	g2 = params(2)
	g3 = params(3)
	s_num = params(4)
	p_num = params(5)
	num = params(6)
	n_i = params(7)
	n_t = params(8)
	row = params(9)
	col = params(10)
	l = params(11)
	
	call boundary(g1, g2, g3)
	
	t1 = 1-g2
	t2 = 2-g2
	h = 1000
	
	allocate(mat(row,col))

	arc1 = 0.5d0*(cal_arc(0d0,1)+cal_arc(pi,1))
	arc2 = 0.5d0*cal_arc(pi*1.5d0,2)
	arc3 = 0.5d0*cal_arc(2d0*pi,3)
	
	do i5 = 1, h
		arc1 = arc1 + cal_arc(0d0+i5*pi/h,1)
		arc2 = arc2 + cal_arc(pi+i5*pi/2d0/h,2)
		arc3 = arc3 + cal_arc(pi*1.5d0+i5*pi/2d0/h,3)
	end do
	
	arc1 = arc1*pi/h
	arc2 = arc2*pi/2d0/h
	arc3 = arc3*pi/2d0/h
	
	arc = arc1 + arc2 + arc3	
	
	Loop_s: do i1 = 1, s_num	
	
		s = 360/s_num*i1
		s = s*pi/180d0
		if (0d0 <= (s*180d0/pi) .AND. (s*180d0/pi) <= 90) then	
			
			r = dsqrt(1d0/(dcos(s)**2+(dsin(s)/g1)**2))
			x = r*dcos(s)
	    	y = r*dsin(s)					
			dy = (x*g1**2)
			dx = -y
			
		else if (91 <= s*180d0/pi .AND. (s*180d0/pi) <= 180) then
			
			r = dsqrt(1d0/(dcos(s)**2+(dsin(s)/g1)**2))
			x = r*dcos(s)
	    	y = r*dsin(s)					
			dy = (x*g1**2)
			dx = -y							
			
		else if (181 <= s*180d0/pi .AND. (s*180d0/pi) <= 270) then
			
			r = dsqrt(1d0/((dcos(s)/(2-g2))**2+(dsin(s)/g3)**2))
			x = r*dcos(s)+1-g2
	    	y = r*dsin(s)		
	    	    	    
	    	dy = (g3**2)*(x-(1-g2))
			dx = -(y*(2-g2)**2)	    	    	
	    	
	    	    
		else
			r = dsqrt(1d0/((dcos(s)/(g2))**2+(dsin(s)/g3)**2))
			x = r*dcos(s)+ 1-g2
	    	y = r*dsin(s)		
	    
	    	dy = (g3**2)*(x-(1-g2))
	    	dx = -(y*(g2)**2)
	    	
	    
		end if	
		
		
		
		x2 = x
		y2 = y
						
		alpha1 = datan2(dy,dx)
		
		Loop_p: do i2 = 1, p_num-1				
			angle = dasin(-1d0+2d0/(p_num)*i2)
            x = x2
            y = y2
        	write(9,11) x,y                                                     
           	alpha = alpha1
			m = alpha + (pi/2) + angle
			m_t = m/pi*180
			m = dtan(m)
			intensity = 100d0
			
			Loop_iteration: do i3 = 1, num            	
			
				IF1: if (dabs(m) > 1e7) then
					if (y >= 0 .and. x <=1-g2) then
						x1 = x					
						y1 = -dsqrt(g3**2*(1-((x-1+g2)**2/(2-g2)**2)))
						bs = 2
						
					else if (y >= 0 .and. x >=1-g2) then
						x1 = x
						y1 = -dsqrt(g3**2*(1-((x-1+g2)**2/(g2)**2)))
						bs = 3

					else if (y < 0) then
						x1 = x
						y1 = dsqrt(g1**2*(1-x**2))
						bs = 1
				
					end if										
					
					write(9,11) x1,y1
				
				else if (dabs(m) < 1e-7) then
					if (y >= 0 .and. x>=0) then
						x1 = -dsqrt(1-(y/g1)**2)
						y1 = y
						bs = 1					
				
					else if (y >=0 .and. x < 0) then
						x1 = dsqrt(1-(y/g1)**2)
						y1 = y
						bs = 1
				
					else if (y < 0 .and. x > 1-g2) then
						x1 = -dsqrt(((2-g2)**2)*(1-(y/g3)**2))+1-g2
						y1 = y
						bs = 2
				
					else if (y < 0 .and. x < 1-g2) then
						x1 = dsqrt(g2**2*(1-(y/g3)**2))+1-g2
						y1 = y
						bs = 3
				
					end if										
					
					write(9,11) x1,y1
			
				
				else
					if (dabs(x) == 1d0) then
						if (x > 0) then
							m1 = y/(x+1)
							m3 = -1e10
						else
							m1 = 1e10
							m3 = -y/(1-x)
						end if
					else																	
						m1 = y/(x+1)
						m3 = -y/(1-x)
					end if
					
					

					if (y >= 0 .AND. x >= 1-g2) then
					
            	        m2 = (y+g3)/(x-1+g2)                    
            	        	
            	        if (m > 0 .AND. m < m1) then
            	        	bs = 1                                    
            	        else if (m > 0 .AND. m1 <= m .AND. m <= m2) then
            	        	bs = 2
            	        else if (m > 0 .AND. m2 < m) then
            	        	bs = 3
            	        else if (m < 0 .AND. dabs(m) <= dabs(m3)) then
            	        	bs = 1
            	        else if (m < 0 .AND. dabs(m) > dabs(m3)) then
            	        	bs = 3
            	        end if
                    
            	    else if (y >= 0 .AND. x < 1-g2) then
                    					
            	        m2 = -(y+g3)/(1-g2-x)
                	    	
            	        if (m > 0 .AND. m >= m1) then
            	        	bs = 2
            	        else if (m > 0 .AND. m < m1) then
            	        	bs = 1
            	        else if (m < 0 .AND. dabs(m) > dabs(m2)) then
            	        	bs = 2
            	        else if (m < 0 .AND. dabs(m) >= dabs(m3) .AND. dabs(m) <= dabs(m2)) then
            	            bs = 3
            	        else if (m < 0 .AND. dabs(m) < dabs(m3)) then
            	        	bs = 1
            	        end if
            	        
            	    else if (y < 0 .AND. x < 1-g2) then
            	        m2 = (-g3-y)/(1-g2-x)
            	        if (m > 0 .AND. m >= m3) then
            	        	bs = 1
            	        else if (m > 0 .AND. m < m3) then
            	            bs = 3
            	        else if (m < 0 .AND. dabs(m) > dabs(m1)) then
            	        	bs = 1
            	        else if (m < 0 .AND. dabs(m) <= dabs(m1) .AND. dabs(m) >= dabs(m2)) then
            	        	bs = 2
            	        else if (m < 0 .AND. dabs(m) < dabs(m2)) then
            	        	bs = 3
            	        end if
            	        
            	    else if (y < 0 .AND. x >= 1-g2) then
            	    	m2 = (y+g3)/(x-1+g2)
            	        if (m > 0 .AND. m > m3) then
            	        	bs = 1
            	        else if (m > 0 .AND. m >= m2 .AND. m <= m3) then
            	        	bs = 3
            	        else if (m > 0 .AND. m < m2) then
            	        	bs = 2
            	        else if (m < 0 .AND. dabs(m) <= dabs(m1)) then
            	        	bs = 2
            	        else if (m < 0 .AND. dabs(m) > dabs(m1)) then
            	        	bs = 1
            	        end if
                    
            	    end if
                    
                	n = y - m*x
                
                	if (bs == 1) then                	                	
                		a = m**2+g1**2
                		b = m*n
                		c = n**2-g1**2
                	                	
                	else if (bs == 2) then                	
                		a = m**2+(g3/t2)**2
                		b = m*n-t1*(g3/t2)**2
                		c = n**2-g3**2+(g3*t1/t2)**2
                	     
                	else if (bs == 3) then
                		a = m**2+(g3/g2)**2
                		b = m*n-t1*(g3/g2)**2
                		c = n**2-g3**2+(g3*t1/g2)**2 
                	                	                	
                	end if                	 
                	
                	e1 = (-b+dsqrt(b**2-a*c))/a
	                e2 = (-b-dsqrt(b**2-a*c))/a
                                       
				    if (dabs(e1) > 1d0) then
						x1 = e2
					else if (dabs(e2) > 1d0) then
						x1 = e1                	                	
                	else if (bs == 2 .AND. e1 > 1-g2) then
                		x1 = e2
                	else if (bs == 2 .AND. e2 > 1-g2) then
                		x1 = e1                		                		
                	else if (bs == 3 .AND. e1 < 1-g2) then
                		x1 = e2
                	else if (bs == 3 .AND. e2 < 1-g2) then
                		x1 = e1
                	else
                		if (dabs(x-e1) < dabs(x-e2)) then
                			x1 = e2
                		else
                			x1 = e1
                		end if                 		                	
                	end if                                										
				
					if (bs == 1) then
						y1 = dsqrt(g1**2*(1-x1**2))
					else if (bs == 2) then
						y1 = -dsqrt(g3**2*(1-((x1-1+g2)/(2-g2))**2))
					else
						y1 = -dsqrt(g3**2*(1-((x1-1+g2)/g2)**2))
					end if									
				
					write(9,11) x1,y1
				
				end if	IF1								
																	
				if (dabs(x1) == 1d0 .OR. dabs(x1) >= 1d0 ) then
					if (x > 0d0) then
						alpha = pi/2d0
					else if (x < 0d0) then
						alpha = pi/2d0*3d0
					end if							
            
				else if (y1 >= 0) then
					dy = (g1**2)*x1
					dx = -dsqrt((g1**2)*(1-x1**2))
					alpha = datan2(dy,dx)												
			
				else if (y1 < 0 .AND. x1 <= 1-g2) then
											    	    	    
		    		dy = (x1-1+g2)*g3**2
					!dx = dsqrt((g3**2)*(1-((x1-1+g2)/t2)**2))
					dx = (2-g2)**2*dsqrt(-(g3**2*(x1+1)*(2*g2+x1-3))/(2-g2)**2)
					alpha = datan2(dy,dx)
	    	    
				else if (y1 < 0 .AND. x1 > 1-g2) then
					dy = (x1-1+g2)*g3**2
					dx = g2**2*dsqrt((g3**2)*(1-((x1-1+g2)/g2)**2))
					alpha = datan2(dy,dx)
				
				else

					print*, "Wrong Decision"

				end if            

			
				k = dtan(alpha + pi/2d0)		
				k1 = (x+x1)/2d0
				k2 = (y+y1)/2d0			
				x_i = (k1/k + k2 + k*x1 - y1)/(1/k+k)
				y_i = k*(x_i-x1)+y1
				v1(1) = x_i-k1
				v1(2) = y_i-k2
				v2(1) = k1-x1
				v2(2) = k2-y1
				v3(1) = 2d0*v1(1) + v2(1)
				v3(2) = 2d0*v1(2) + v2(2)
						
				if (dabs(k) > 1e10) then				
					if (k > 0) then
						k = 1e10
					else
						k = -1e10
					end if
					x_i = (k1/k + k2 + k*x1 - y1)/(1/k+k)
					y_i = k*(x_i-x1)+y1
					v1(1) = x_i-k1
					v1(2) = y_i-k2
					v2(1) = k1-x1
					v2(2) = k2-y1
					v3(1) = 2d0*v1(1) + v2(1)
					v3(2) = 2d0*v1(2) + v2(2)
					m = v3(2)/v3(1)	

				else if (dabs(k) < 1e-10) then 				
					if (k > 0) then					
						k = 1e-10
					else 					
						k = -1e-10		
					end if
					x_i = (k1/k + k2 + k*x1 - y1)/(1/k+k)
					y_i = k*(x_i-x1)+y1
					v1(1) = x_i-k1
					v1(2) = y_i-k2
					v2(1) = k1-x1
					v2(2) = k2-y1
					v3(1) = 2d0*v1(1) + v2(1)
					v3(2) = 2d0*v1(2) + v2(2)
					m = v3(2)/v3(1)					
			
				else
					m = v3(2)/v3(1)	
				end if									
			
				v4(1) = x_i - x1
				v4(2) = y_i - y1
				p_psos = v4(1)*v2(2)-v4(2)*v2(1)

				if (v2(1) == 0d0 .AND. v2(2) == 0d0) then
					p_psos = 0d0
				else if (v4(1) == 0d0 .AND. v4(2) == 0d0) then
					p_psos = 0d0
				else
					p_psos = p_psos/(dsqrt(v4(1)**2+v4(2)**2)*dsqrt(v2(1)**2+v2(2)**2))
				end if

				if (bs == 1) then

					theta_arc = dacos(x1)

					s_psos = 0.5d0*(cal_arc(0d0,1)+cal_arc(theta_arc,1))
					do i5 = 1, h
						s_psos = s_psos + cal_arc(0d0+i5*theta_arc/h,1)
					end do
					s_psos = s_psos*theta_arc/h

				else if (bs == 2) then

					theta_arc = dacos(x1)

					theta_arc = 2d0*pi - dacos(x1)
					s_psos = 0.5d0*(cal_arc(pi,2)+cal_arc(theta_arc,2))
					do i5 = 1, h
						s_psos = s_psos + cal_arc(pi+i5*theta_arc/h,2)
					end do				
					s_psos = s_psos*(theta_arc-pi)/h
					s_psos = s_psos + arc1
				else 

					theta_arc = 2d0*pi - dacos(x1)				 				
					s_psos = 0.5d0*(cal_arc(1.5d0*pi,3)+cal_arc(theta_arc,3))
					do i5 = 1, h
						s_psos = s_psos + cal_arc(1.5d0*pi+i5*theta_arc/h,3)
					end do								
					s_psos = s_psos*(theta_arc-pi*1.5d0)/h
					s_psos = s_psos + arc1 + arc2

				end if

				s_psos = s_psos/arc

				theta_i = dasin(dabs(p_psos))

				if (theta_i <= n_t/n_i) then
					
					!SPD
					theta_t = dasin(n_i/n_t*dabs(p_psos))
					r_s = (n_i*dcos(theta_i)-n_t*dcos(theta_t))/(n_i*dcos(theta_i)+n_t*dcos(theta_t))				
					ref = r_s**2
					emi = 1d0-ref
					intensity = intensity*ref

					!Far-field
					if (i3 >= 0.9d0*num) then
					v_f(1) = -v4(1)
					v_f(2) = -v4(2)

					v_f(1) = v_f(1) + v1(1)/dsqrt(v1(1)**2+v1(2)**2)*dtan(theta_t)*dsqrt(v4(1)**2+v4(2)**2)
					v_f(2) = v_f(2) + v1(2)/dsqrt(v1(1)**2+v1(2)**2)*dtan(theta_t)*dsqrt(v4(1)**2+v4(2)**2)

					m_f = v_f(2)/v_f(1)

					q = y1 - m_f*x1

					a_f = m_f**2 + 1d0
					b_f = m_f*q
					c_f = q**2 - l**2

					e_f1 = (-b_f+dsqrt(b_f**2-a_f*c_f))/a_f
	        	   	e_f2 = (-b_f-dsqrt(b_f**2-a_f*c_f))/a_f

					if (v_f(1) >= 0) then
						x_f = e_f1
					else
						x_f = e_f2
					end if

					y_f = m_f*x_f + q
				
					angle_f = dacos(x_f/l)

					if (y_f < 0) then
						angle_f = 2*pi - angle_f	
					end if				

					angle_f = angle_f/pi*180d0

					angle_f = floor(angle_f) + 1

					mat2(angle_f) = mat2(angle_f) + intensity

					end if
				end if 
				

				if (i3 >= 0.9d0*num) then
				
					tmp_s = s_psos
					tmp_p = p_psos

					do i = 1, col
						tmp_s = tmp_s - 1d0/col
						if (tmp_s < 0) then
							tmp_r = i
							EXIT
						end if
					end do		

					do j = 1, row
						tmp_p = tmp_p + 2d0/row
						if (tmp_p > 1) then
							tmp_c = j
							EXIT
						end if
					end do 

				 
					mat(tmp_r, tmp_c) = mat(tmp_r, tmp_c) + intensity
				
				end if
				
				

				write(1111,11) s_psos, p_psos
				write(321,11) x1, y1
				write(321,11) x_f, y_f
				write(321,11)
				
			
				x = x1
				y = y1
			
        	end do Loop_iteration											


		end do Loop_p

		print*, "Progress:", i1*i2*i3/(s_num*p_num*num)*100, "%"
	
	end do Loop_s
	
	do j = 0, row-1
		tmp_y = 1d0 - 1d0/row - j*2d0/row

		do i = 0, col-1
			tmp_x = 1d0/2d0/col + i*1d0/col
			write(4321,11) tmp_x, tmp_y, mat(i+1, j+1)
		end do
		write(4321,11) 

	end do

	do i = 1, 361
		angle_w = i
		write(123,11) angle_w, mat2(i) 
	end do


	close(9873)
	close(9)
	close(1111)
	close(1234)
	close(4321)
	close(123)
	close(321)
	deallocate(mat)
				
	contains
												
	subroutine boundary(g1, g2, g3)
       implicit none
       real(8), INTENT(IN) :: g1, g2, g3
       real(8) :: x, y, r, pi, theta
       integer :: i0
       
       open(unit = 10,file = "./boundary.dat")
       11 format(X,4E)
                     
       pi = 2d0*dasin(1d0)
       
       
       do i0 = 1,180
       theta = real(i0,8)*pi/180d0
       r = dsqrt(1d0/(dcos(theta)**2+(dsin(theta)/g1)**2))
       x = r*dcos(theta)
       y = r*dsin(theta)
       write(10,11) x,y
       end do
       
       do i0 = 181,270
       theta = real(i0,8)*pi/180d0
       r = dsqrt(1d0/((dcos(theta)/(2-g2))**2+(dsin(theta)/g3)**2))
       x = r*dcos(theta)
       y = r*dsin(theta)
       write(10,11) x+(1-g2),y
       end do
       
       do i0 = 271,360
       theta = real(i0,8)*pi/180d0
       r = dsqrt(1d0/((dcos(theta)/(g2))**2+(dsin(theta)/g3)**2))
       x = r*dcos(theta)
       y = r*dsin(theta)
       write(10,11) x+(1-g2),y
       end do
       
       close(10)
       
       end subroutine boundary 
		
	function radius(theta, bs, diff)
	
		implicit none
		
		real(8), INTENT(IN) :: theta 
		integer, INTENT(IN) :: bs, diff
		real(8) :: radius, g1, g2, g3, e 
		common /coeff/ g1,g2,g3
		
		if (bs == 1) then
			a = 1d0
			b = g1
		else if (bs == 2) then
			a = 2d0-g2
			b = g3
		else 
			a = g2
			b = g3
		end if
		if (a >= b) then
			e = dsqrt(1d0-(b/a)**2)
			if (diff == 0) then
				radius = b/dsqrt(1d0-(e*dcos(theta))**2)
			else
				radius = -e**2*b*dsin(theta)*dcos(theta)/(1d0-(e*dcos(theta))**2)**(3d0/2d0)
			end if
		else
			e = dsqrt(1d0-(a/b)**2)
			if (diff == 0) then
				radius = a/dsqrt(1d0-(e*dcos(theta))**2)
			else
				radius = -e**2*a*dsin(theta)*dcos(theta)/(1d0-(e*dcos(theta))**2)**(3d0/2d0)
			end if
		end if
		
		
		
	end function radius
	
	function cal_arc(theta, bs)
	
		implicit none
		
		real(8), INTENT(IN) :: theta
		integer :: bs
		real(8) :: cal_arc, dxdt, dydt
		
		dxdt = radius(theta, bs, 1)*dcos(theta) - radius(theta, bs, 0)*dsin(theta)
		dydt = radius(theta, bs, 1)*dsin(theta) + radius(theta, bs, 0)*dcos(theta)
		
		cal_arc = dsqrt(dxdt**2+dydt**2)								
	
	end function cal_arc
			


	end program nrm
	
	

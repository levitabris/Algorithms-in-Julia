```
  This is a direct Julia port of Mark G. Johnson's C implementation of the Hooke and Jeeves algorithm based on Sean R. Johnson's Python version on July 7, 2013. The original source code can be found at https://searchcode.com/file/66639051/hooke_jeeves_bounded.py.
  
  ===== Original Copyright Information ======
  
  Nonlinear Optimization using the algorithm of Hooke and Jeeves  
  12 February 1994	author: Mark G. Johnson 

  The author of this software is M.G. Johnson.		   
  Permission to use, copy, modify, and distribute this software
  for any purpose without fee is hereby granted, provided that
  this entire notice is included in all copies of any software
  which is or includes a copy or modification of this software
  and in all copies of the supporting documentation for such
  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT
  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE
  AUTHOR NOR AT&T MAKE ANY REPRESENTATION OR WARRANTY OF ANY
  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS
  FITNESS FOR ANY PARTICULAR PURPOSE.
 
  ==========================================
```



function rosenbrock(x)
    a = x[1]
    b = x[2]
    ((1.0 - a) ^ 2) + (100.0 * (b - (a ^ 2)) ^ 2)
end

function _hooke_best_nearby(f, delta, point, prevbest, bounds=nothing, args=[])
    z = copy(point)
    minf = prevbest
    ftmp = 0.0
    fev = 0
    
    for i = 1:length(points)
        #see if moving point in the positive delta direction decreases the 
        z[i] = _value_in_bounds(point[i] + delta[i], bounds[i][1], bounds[i][2])
        
        ftmp = f(z, args...)
        fev += 1
        if ftmp < minf
            minf = ftmp
        else
            #if not, try moving it in the other direction
            z[i] = _value_in_bounds(point[i] - delta[i], bounds[i][1], bounds[i][2])
            
            ftmp = f(z, args...)
            fev += 1
            if ftmp < minf
                minf = ftmp
            else
                #if moving the point in both delta directions result in no improvement, then just keep the point where it is
                z[i] = point[i]
            end
        end
    end
    for i = 1:length(x)
        point[i] = z[i]
        return (minf, fev)
    end
end

function _value_in_bounds(val, low, high)
    if val < low
        return low
    elseif val > high
        return high
    else
        return val
    end
end

function _point_in_bounds(point, bounds)
    for i = 1:(size(point)[1])
        if point[i] < bounds[i][1]
            point[i] = bounds[i][1]
        elseif point[i] > bounds[i][2]
            point[i] = bounds[i][2]
        end
    end
end

function _is_point_in_bounds(point, bounds)
    out = true
    for i = 1:length(point)
        if point[i] < bounds[i][1]
            out = false
        elseif point[i] > bounds[i][2]
            out = false
        end
    end
    return out
end

function hooke(f, startpt, bounds=nothing, rho=0.5, epsilon=1e-6, itermax=5000, args=[])

    result = Dict(
        "success"=> true,
        "message" => "Success"
    )
    
    delta = zeros(length(startpt))
    endpt = zeros(length(startpt))

    if isnothing(bounds)
        # if bounds is none, make it none for all (it will be converted to below)
        bounds = [[nothing, nothing] for x in startpt]
    else
        @bp
        bounds = [[x[1], x[2]] for x in bounds] #make it so it wont update the original
    end

    startpt = copy(startpt) #make it so it wont update the original
    
    fmin = nothing
    nfev = 0
    iters = 0
    
    for bound in bounds
        if isnothing(bound[1])
            bound[1] = Inf
        else
            bound[1] = float(bound[1])
        end
        if isnothing(bound[2])
            bound[2] = Inf
        else
            bound[2] = float(bound[2])
        end
    end
    try
        # shift 
        _point_in_bounds(startpt, bounds) #shift startpt so it is within the bounds
        
        xbefore = copy(startpt)
        newx = copy(startpt)
        for i = 1:length(startpt)
            delta[i] = abs(startpt[i] * rho)
            if (delta[i] == 0.0)
                # we always want a non-zero delta because otherwise we'd just be checking the same point over and over
                # and wouldn't find a minimum
                delta[i] = rho
            end
        end

        steplength = rho

        fbefore = f(newx, args...)
        nfev += 1
        
        newf = fbefore
        fmin = newf
        while ((iters < itermax) & (steplength > epsilon))
            iters += 1
            # print "after %5d , f(x) = %.4le at" % (funevals, fbefore)
            # for j in range(len(startpt)):
                #print "   x[%2d] = %4le" % (j, xbefore[j])
            # pass
            
            ##/* find best new point, one coord at a time */
            newx = copy(xbefore)
            (newf, evals) = _hooke_best_nearby(f, delta, newx, fbefore, bounds, args)
            
            nfev += evals
            ##/* if we made some improvements, pursue that direction */
            keep = 1
            while (newf < fbefore) & (keep == 1)
                fmin = newf
                for i = 1:length(startpt)
                    ##/* firstly, arrange the sign of delta[] */
                    if newx[i] <= xbefore[i]
                        delta[i] = -abs(delta[i])
                    else
                        delta[i] = abs(delta[i])
                    end
                    ## /* now, move further in this direction */
                    tmp = xbefore[i]
                    xbefore[i] = newx[i]
                    newx[i] = _value_in_bounds(newx[i] + newx[i] - tmp, bounds[i][1], bounds[i][2])
                end
                fbefore = newf
                (newf, evals) = _hooke_best_nearby(f, delta, newx, fbefore, bounds, args)
                nfev += evals
                ##/* if the further (optimistic) move was bad.... */
                if (newf >= fbefore)
                    break
                end
                ## /* make sure that the differences between the new */
                ## /* and the old points are due to actual */
                ## /* displacements; beware of roundoff errors that */
                ## /* might cause newf < fbefore */
                keep = 0
                for i = 1:length(startpt)
                    keep = 1
                    if ( abs(newx[i] - xbefore[i]) > (0.5 * abs(delta[i])) )
                        break
                    else
                        keep = 0
                    end
                end
            end
            if ((steplength >= epsilon) & (newf >= fbefore))
                steplength = steplength * rho
                delta = [x * rho for x in delta]
            end
        end
        for x = 1:length(xbefore)
            endpt[x] = xbefore[x]
        end
    catch e
        warn("Exception: ", e)
        result["success"] = false
        result["message"] = e
    finally
        result["nit"] = iters
        result["fevals"] = nfev
        result["fun"] = fmin
        result["x"] = endpt
        return result
    end
end

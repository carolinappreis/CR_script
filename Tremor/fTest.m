        % Test if non-linear part of the model function is zero. Do one of
        % the following:
        % (1) If full model has an intercept and number of estimated coefficients
        %     is > 1 then test full model against the intercept only model.  
        % (2) Otherwise, test full model against a zero model.
        function [f,p,emptyNullModel,hasIntercept] = fTest(model)
        % Variable emptyNullModel is:
        %  -  true if our null model is the zero model. 
        %  - false if our null model is the intercept only model.            
            
            % (1) Inspect unweighted Jacobian and figure out if there is an
            % intercept (i.e., if if Junw_r has a constant column).
            [~,Junw_r] = create_J_r(model);
            Jmin = min(Junw_r,[],1);
            Jmax = max(Junw_r,[],1);
            hasIntercept = any( abs(Jmax-Jmin) <= sqrt(eps(class(Junw_r))) * (abs(Jmax) + abs(Jmin)) ); 
            if hasIntercept &&  (model.NumEstimatedCoefficients > 1)
                % (2) Compare full model vs. intercept only model.
                emptyNullModel = false;
                nobs = model.NumObservations;
                ssr = max(model.SST - model.SSE,0);
                dfr = model.NumEstimatedCoefficients - 1;
                dfe = nobs - 1 - dfr;
                f = (ssr./dfr) / (model.SSE/dfe);
                p = fcdf(1./f,dfe,dfr); % upper tail
            else
                % (2) Compare full model vs. zero model.
            	emptyNullModel = true;
                ssr = max(model.SST0 - model.SSE,0);
                dfr = model.NumEstimatedCoefficients;
                dfe = model.NumObservations - model.NumEstimatedCoefficients;
                f = (ssr./dfr) / (model.SSE/dfe);
                p = fcdf(1./f,dfe,dfr); % upper tail
            end                        
        end 
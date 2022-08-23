function [SM, EN, SYr, SXr, SE, SG, OPP] = Extrema_Single_Symm_Scan(img)
% "Search of Single Extrema using Symmetric scan-line order", made by Tuan Nguyen, Tsvetkov, 2019
%
% img: input graysacle or RGB image
% SE: stack of number of local extrema
% SG: stack of extreme values
% BH: check matrix of maxima
% BL: check matrix of minima
% OPP: Number of comparisons per pixel
% SM: Label matrix of local extrema
% EN: number of extrema
% SYr, SXr: coordinates, seed values of first pixels of extrema.
% Memory Use (bytes): ~72*N, N-Number of pixels
% Tuan Q. Pham. Non-maximum Suppression Using fewer than 2 Comparisons per Pixel. Advanced Concepts for Intelligent Vision
% Systems (ACIVS) 12, Sydney, Australia, December 13-16 (2010), pp. 438–451.
%
% Convert RGB to Grayscale image
if ndims(img)==3
    img = rgb2gray(img);
end
img=double(img);
% Input paramaters
[Y,X] = size(img);
SE=zeros(Y*X,2);
SYr=zeros(1,Y*X); SXr=zeros(1,Y*X); ENr=0;
N = 0;
EN = 0;
SM = zeros(Y,X);
for y =2:Y-1
    for x=2:X-1
        N=N+10;
        % Scanning Mask: [6 7 8; 2 0 1; 5 4 3]
        S6=img(y-1,x-1); S7=img(y-1,x); S8=img(y-1,x+1);
        S2=img(y,x-1);   S0=img(y,x);   S1=img(y,x+1);
        S5=img(y+1,x-1); S4=img(y+1,x); S3=img(y+1,x+1);
            %1.Find local minima using non-maximum suppression
                if S1>=S0
                    N=N+1;
                    if S1>S0
                        N=N+1;
                        if S2<=S0
                            N=N+1;
                        else%S2>S0
                            N=N+1;
                            if S3<=S0
                                N=N+1;
                            else%S3>S0
                                N=N+1;
                                if S4<=S0
                                    N=N+1;
                                else%S4>S0
                                    N=N+1;
                                    if S5<=S0
                                        N=N+1;
                                    else%S5>S0
                                        N=N+1;
                                        if S6<=S0
                                            N=N+1;
                                        else%S6>S0
                                            N=N+1;
                                            if S7<=S0
                                                N=N+1;
                                            else%S7>S0
                                                N=N+1;
                                                if S8<=S0
                                                    N=N+1;
                                                else%S8>S0
                                                    N=N+6;
                                                    EN=EN+1;
                                                    SM(y,x)=-EN;
                                                    ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
                                                    SE(EN,1)=1;
                                                    SE(EN,2)=img(y,x);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                else%S1<S0
                    N=N+1;
                    if S2>=S0
                        N=N+1;
                    else%S2<S0
                        N=N+1;
                        if S3>=S0
                            N=N+1;
                        else%S3<S0
                            N=N+1;
                            if S4>=S0
                                N=N+1;
                            else%S4<S0
                                N=N+1;
                                if S5>=S0
                                    N=N+1;
                                else%S5<S0
                                    N=N+1;
                                    if S6>=S0
                                       N=N+1;
                                    else%S6<S0
                                        N=N+1;
                                        if S7>=S0
                                            N=N+1;
                                        else%S7<S0
                                            N=N+1;
                                            if S8>=S0
                                                N=N+1;
                                            else%S8<S0
                                                N=N+6;
                                                EN=EN+1;
                                                SM(y,x)=EN;
                                                ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
                                                SE(EN,1)=1;
                                                SE(EN,2)=img(y,x);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
          %
    end
end
%Result
%E=sign(EM);
% Compute the comparisons per pixel
OPP = N/(Y*X);
SG=round(SE);
end
%%
% function [SM, EN, SYr, SXr, SE, SG, OPP] = Extrema_Single_Symm_Scan(img)
% % "Search of Single Extrema using Symmetric scan-line order", made by Tuan Nguyen, Tsvetkov, 2019
% %
% % img: input graysacle or RGB image
% % SE: stack of number of local extrema
% % SG: stack of extreme values
% % BH: check matrix of maxima
% % BL: check matrix of minima
% % OPP: Number of comparisons per pixel
% % SM: Label matrix of local extrema
% % EN: number of extrema
% % SYr, SXr: coordinates, seed values of first pixels of extrema.
% % Memory Use (bytes): ~72*N, N-Number of pixels
% % Tuan Q. Pham. Non-maximum Suppression Using fewer than 2 Comparisons per Pixel. Advanced Concepts for Intelligent Vision
% % Systems (ACIVS) 12, Sydney, Australia, December 13-16 (2010), pp. 438–451.
% %
% % Convert RGB to Grayscale image
% if ndims(img)==3
%     img = rgb2gray(img);
% end
% img=double(img);
% % Input paramaters
% [Y,X] = size(img);
% SE=zeros(Y*X,2);
% SYr=zeros(1,Y*X); SXr=zeros(1,Y*X); ENr=0;
% BH = zeros(Y,X);
% BL = zeros(Y,X);
% N = 0;
% EN = 0;
% SM = zeros(Y,X);
% for y =2:Y-1
%     for x=2:X-1
%         N=N+10;
%         % Scanning Mask: [6 7 8; 2 0 1; 5 4 3]
%         S6=img(y-1,x-1); S7=img(y-1,x); S8=img(y-1,x+1);
%         S2=img(y,x-1);   S0=img(y,x);   S1=img(y,x+1);
%         S5=img(y+1,x-1); S4=img(y+1,x); S3=img(y+1,x+1);
%             %1.Find local minima using non-maximum suppression
%             BH(y,x) = 1;
%                 if S1>=S0
%                     BH(y,x)=0; N=N+2;
%                     if S1>S0
%                         BL(y,x)=1; N=N+2;
%                         if S2<=S0
%                             BL(y,x)=0; N=N+2;
%                         else%S2>S0
%                             N=N+1;
%                             if S3<=S0
%                                 BL(y,x)=0; N=N+2;
%                             else%S3>S0
%                                 N=N+1;
%                                 if S4<=S0
%                                     BL(y,x)=0;N=N+2;
%                                 else%S4>S0
%                                     N=N+1;
%                                     if S5<=S0
%                                         BL(y,x)=0;N=N+2;
%                                     else%S5>S0
%                                         N=N+1;
%                                         if S6<=S0
%                                             BL(y,x)=0;N=N+2;
%                                         else%S6>S0
%                                             N=N+1;
%                                             if S7<=S0
%                                                 BL(y,x)=0;N=N+2;
%                                             else%S7>S0
%                                                 N=N+1;
%                                                 if S8<=S0
%                                                     BL(y,x)=0;N=N+2;
%                                                 else%S8>S0
%                                                     N=N+6;
%                                                     EN=EN+1;
%                                                     SM(y,x)=EN;
%                                                     ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
%                                                     SE(EN,1)=1;
%                                                     SE(EN,2)=img(y,x);
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 else%S1<S0
%                     N=N+1;
%                     if S2>=S0
%                         BH(y,x)=0;N=N+2;
%                     else%S2<S0
%                         N=N+1;
%                         if S3>=S0
%                             BH(y,x)=0;N=N+2;
%                         else%S3<S0
%                             N=N+1;
%                             if S4>=S0
%                                 BH(y,x)=0;N=N+2;
%                             else%S4<S0
%                                 N=N+1;
%                                 if S5>=S0
%                                     BH(y,x)=0;N=N+2;
%                                 else%S5<S0
%                                     N=N+1;
%                                     if S6>=S0
%                                         BH(y,x)=0;N=N+2;
%                                     else%S6<S0
%                                         N=N+1;
%                                         if S7>=S0
%                                             BH(y,x)=0;N=N+2;
%                                         else%S7<S0
%                                             N=N+1;
%                                             if S8>=S0
%                                                 BH(y,x)=0;N=N+2;
%                                             else%S8<S0
%                                                 N=N+6;
%                                                 EN=EN+1;
%                                                 SM(y,x)=-EN;
%                                                 ENr=ENr+1; SYr(ENr)=y; SXr(ENr)=x;
%                                                 SE(EN,1)=1;
%                                                 SE(EN,2)=img(y,x);
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%           %
%     end
% end
% %Result
% %E=sign(EM);
% % Compute the comparisons per pixel
% OPP = N/(Y*X);
% SG=round(SE);
% end
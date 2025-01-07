

function isOpen = isExcelFileOpen(filename,skipFlag)
% isOpen = isExcelFileOpen(filename,skipFlag)
%
% filename = name of the excel file...
% skipFlag: if true, this won't delay until the sheet is closed, it will
% just report that it is open
% 
% tests that the excel file is writable before trying it with the main code
% can skip the file address, but it is faster to give the address 
% 
% ES Dickinson


  if nargin==0
      [~, ~, filename] = load_QuadBowlExperiments;
  end
    
  if nargin<2
      skipFlag = false;
  end

  isOpen = true;

  while isOpen
      try
        writecell({'test'},filename,'Sheet','Sheet1','Range', 'A1')
        isOpen = false;
      catch
            if ~skipFlag
                h = warndlg('Close the Excel file and try again', 'Close Excel File');
                uiwait(h)
            else 
                return
            end
      end
    
  end



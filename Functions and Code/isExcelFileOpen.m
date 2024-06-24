

function isOpen = isExcelFileOpen(filename)
% isOpen = isExcelFileOpen(filename)
%
% tests that the excel file is writable before trying it with the main code
% can skip the file address, but it is faster to give the address 
%
  if nargin==0
      [~, ~, filename] = load_QuadBowlExperiments;
  end

  isOpen = true;

  while isOpen
      try
        writecell({'test'},filename,'Sheet','Sheet1','Range', 'A1')
        isOpen = false;
      catch
            h = warndlg('Close the Excel file and try again', 'Close Excel File');
            uiwait(h)
      end
    
  end


